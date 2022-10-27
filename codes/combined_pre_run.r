#--------------Setup----------------------

#loading required packages (You should install all of them. Some need to be installed by Bioconductor)

library(igraph)
library(doParallel)
library(parallel)
library(data.table)
library(pbapply)
library(stringr)
library(reshape2)
library(iterators)
library(cluster)
library(limma)
library(foreign)
library(MASS)
library(topGO)
library(colorRamps)
library(GO.db)
library(ROCR)
library(survival)
library(robustbase)
library(ConsensusClusterPlus)
library(RColorBrewer)

#data.frame will no longer convert string to factor
options(stringsAsFactors = F)

#------------------Predefined Functions------------------------------

#Export table with row names and column names
export.1=function(t,filename=NULL){
  if(is.null(filename)) filename=paste0(deparse(substitute(t)),'.txt')
  write.table(t,filename,quote=F,sep='\t',na="None")
}

#Export table without rownames and colnames
export.2=function(t,filename=NULL){
  if(is.null(filename)) filename=paste0(deparse(substitute(t)),'.txt')
  write.table(t,filename,quote=F,sep='\t',na="None",row.names=F,col.names=F)
}

#Import table with row names and column names
#s: file directory
#r: assign (return) the table to a variable or create a new variable based on the file name
#nas: marker for "Not Available"
#sep: delimiter
import.1=function(s,r=F,nas="",sep="\t"){
  if(r==T){
    return(read.delim(s,na.strings=nas,quote="",sep=sep,stringsAsFactors=F))
  }else{
    s1=strsplit(s,'/',T)[[1]]
    s1=s1[length(s1)]
    assign(substring(s1,1,nchar(s1)-4),read.delim(s,na.strings=nas,quote="",stringsAsFactors=F,sep=sep),envir=.GlobalEnv)
  }
}

#Import table without rownames and colnames
import.2=function(s,r=F,nas="",sep='\t'){
  if(r==T){
    return(read.delim(s,header=F,na.strings=nas,sep=sep,quote="",stringsAsFactors=F))
  }else{
    s1=strsplit(s,'/',T)[[1]]
    s1=s1[length(s1)]
    assign(substring(s1,1,nchar(s1)-4),read.delim(s,header=F,na.strings=nas,sep=sep,quote="",stringsAsFactors=F),envir=.GlobalEnv)
  }
}

#Table replacement function
#t: a data.frame that needs be replaced
#d: the dictionary with the first column as key and the second column as the value
#cn: the colnames or indices of the columns need to be replaced
#war: show warning if there are duplicated key in d
replace.1=function(t,d,cn=NULL,war=T){
  t=data.table(t)
  setkey(t,NULL)
  d=data.table(d)
  d=d[,1:2,with=F]
  setnames(d,1:2,c('rep.key','rep.value'))
  cn.o=copy(names(t))
  if(is.numeric(cn)) cn=cn.o[cn] else if(is.null(cn)) cn=cn.o
  if(is.factor(d[[1]])) d[,names(d)[1]:=levels(d[[1]])[as.integer(d[[1]])]]
  key.class=class(d[[1]])
  t[, cn:=lapply(.SD,function(x){
    if(is.factor(x)) as(levels(x)[as.integer(x)],key.class) else as(x,key.class)}),
    .SDcols=cn,with=F]
  setkeyv(d,names(d)[1])
  if(any(duplicated(d))){
    if(war) warning("key is not unique",call.=F)
    for (i in cn) {
      cn.n=c(i,setdiff(cn.o,i))
      setcolorder(t,cn.n)
      estr=paste0("list(",i,"=",colnames(d)[2],if(length(cn.n)==1) "" else ",",paste0(cn.n[-1],collapse=","),")")
      estr=parse(text=estr)
      t=d[t,eval(estr),allow.cartesian=T]
    }
    setcolorder(t,cn.o)
  }else{
    t[, cn:=lapply(.SD,function(x) d[J(x)][[names(d)[2]]]),.SDcols=cn,with=F]
  }
  return(t)
}

#Test whether a subnetwork is connected
#g: Network
#vs: Subnetwork nodes
test.connected=function(g,vs){
  g=induced.subgraph(g,vs)
  return(is.connected(g))
}

#Test whether a subnetwork is connected by adjacency matrix
#adj.m: Network adjacency matrix
#x: Subnetwork nodes
#n: Subnetwork size (length(x))
test.connected2=function(adj.m,x,n){
  m1=adj.m[x,x]
  temp=colSums(m1)
  l1=diag(temp)-m1
  return(round(eigen(l1,symmetric=T,only.values=T)$values[n-1],5)>0)
}

#Vector standardization
getZ2=function(x){
  m=mean(x)
  s=sd(x)
  if(s==0){
    return(rep(0,length(x)))
  }else{
    return((x-m)/s)
  }
}

#Calculate euclidean distance between two vectors
euclidean.dist=function(x1,x2){
  return(sqrt(sum((x1-x2)^2)))
}

#Calculate absolute pearson distance between two vectors
abspearson.dist=function(x1,x2){
  return(1-abs(sum(x1*x2)/sqrt(sum(x1^2)*sum(x2^2))))
}

#Remove islands in a network
#x: Network
#min.v: Maximum size of an island
remove.graph.islands=function(x,min.v=3){
  x1=decompose.graph(x,min.vertices = min.v)
  valid.vs=unlist(lapply(x1,function(y)V(y)$name))
  if(length(x1)==0) return(NULL) else return(induced.subgraph(x,valid.vs))
}

#Define color scales
heat.col=colorRampPalette(c("blue","white","red"),space='rgb')
heat.col3=heat.col(299)

#Remove outliers. See adjboxStats package for detail
remove.outliers=function(s.m,coef1=1.5,a1=-4,b1=3){
  sm2=s.m
  sm3=as.numeric(sm2)
  qs=adjboxStats(sm3,coef=coef1,a=a1,b=b1)$stats[c(1,5)]
  sm2[sm2<qs[1]]=qs[1]
  sm2[sm2>qs[2]]=qs[2]
  return(sm2)
}

#Generate Heatmap (function: heatmap.31)
source('heatmap31.R')


#k-fold cross-validation sampling
createFolds=function (y, k = 10, list = TRUE, returnTrain = FALSE)
{
    if (is.numeric(y)) {
        cuts <- floor(length(y)/k)
        if (cuts < 2)
            cuts <- 2
        if (cuts > 5)
            cuts <- 5
        breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
        y <- cut(y, breaks, include.lowest = TRUE)
    }
    if (k < length(y)) {
        y <- factor(as.character(y))
        numInClass <- table(y)
        foldVector <- vector(mode = "integer", length(y))
        for (i in 1:length(numInClass)) {
            min_reps <- numInClass[i]%/%k
            if (min_reps > 0) {
                spares <- numInClass[i]%%k
                seqVector <- rep(1:k, min_reps)
                if (spares > 0)
                  seqVector <- c(seqVector, sample(1:k, spares))
                foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
            }
            else {
                foldVector[which(y == names(numInClass)[i])] <- sample(1:k,
                  size = numInClass[i])
            }
        }
    }
    else foldVector <- seq(along = y)
    if (list) {
        out <- split(seq(along = y), foldVector)
        names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))),
            sep = "")
        if (returnTrain)
            out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
    }
    else out <- foldVector
    out
}

#Generate color vector based on labels and color scale
make.color=function(labels,color.fun=rainbow){
  un.elements=sort(unique(labels),na.last=T)
  if(any(is.na(un.elements))){
    col.length=length(un.elements)-1
    cols=c(color.fun(col.length),'#BEBEBE')
    legend=data.frame(Label=c(un.elements[-length(un.elements)],'N/A'),
                      Color=cols,stringsAsFactors=F)
    t1=labels
    t1[is.na(t1)]='--N/A--'
    color=t1
    Freq=c()
    for(i in 1:(length(un.elements)-1)){
      c.l=un.elements[-length(un.elements)][i]
      Freq=c(Freq,sum(color==c.l))
      color[color==c.l]=legend[legend$Label==c.l,]$Color
    }
    Freq=c(Freq,sum(color=='--N/A--'))
    color[color=='--N/A--']=legend[legend$Label=='N/A',]$Color
    legend$Freq=Freq
    return(list(legend=legend,color=color))
  }else{
    col.length=length(un.elements)
    cols=color.fun(col.length)
    legend=data.frame(Label=un.elements,
                      Color=cols,stringsAsFactors=F)
    color=labels
    Freq=c()
    for(i in 1:length(un.elements)){
      c.l=un.elements[i]
      Freq=c(Freq,sum(color==c.l))
      color[color==c.l]=legend[legend$Label==c.l,]$Color
    }
    legend$Freq=Freq
    return(list(legend=legend,color=color))
  }
}

#Plot legend
plot.legend=function(legend,...){
  legend=legend$legend
  par(mar=c(5.1,10,4.1,2.1))
  barplot(legend$Freq,names.arg=legend$Label,col=legend$Color,horiz=T,
          las=2,xlab='Frequency',...)
  par(mar=c(5.1,4.1,4.1,2.1))
}

#Generate color based on expression values and color scale
get.color2=function(gs,samples,exp.m=Sbn.Exp,max.lim=1.5,color.f=heat.col){
  ms=sapply(gs,function(x)mean(exp.m[x,samples]))
  cs=color.f(299)
  ms[ms>max.lim]=max.lim
  ms[ms<(-max.lim)]=-max.lim
  return(cs[round((ms+max.lim)/2/max.lim*298)+1])
}

#Generate custom distance matrix
#my.matrix: Matrix with columns as observations
#my.function: distance function
custom.corm <- function(my.matrix, my.function) {
  n <- NCOL(my.matrix)
  mat <- matrix(0, ncol = n, nrow = n)
  colnames(mat) <- rownames(mat) <- colnames(my.matrix)
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      mat[i,j] <- my.function(my.matrix[,i],my.matrix[,j])
    }}
  return(mat)
}

#------------------GSEA functions---------------------------------------
#Setup data for gene enrichment analysis
#gs: Significant gene list
#allgs: All gene list
#ns: Minimum number of genes allowed in a term
#mapping.base: GO dictionary
newGOdata=function(gs,allgs,ns=3,mapping.base='org.Hs.eg.db'){
  ags=rep(0,length(allgs))
  names(ags)=allgs
  ags[allgs %in% gs]=1
  data=new("topGOdata",description='Subnetwork Genes',ontology="BP",allGenes=ags,geneSel=function(x)x==1,nodeSize=ns,annot=annFUN.org,mapping=mapping.base,ID="symbol")
  return(data)
}

#Gene enrichment analysis
#data: data generated by newGOdata
#alg: Algorithm
#test: Significance test
GOtest=function(data,alg="parentchild",test="fisher"){
  result=runTest(data, algorithm=alg, statistic = test)
  return(result)
}

getSigGoTerms=function(gs,allgs,ns=3,alg="lea",test="fisher",mapping.base='org.Hs.eg.db'){
  gdata=newGOdata(gs,allgs,ns,mapping.base)
  result=GOtest(gdata,alg,test)
  r.table=GenTable(gdata,Fisher=result,topNodes=20,numChar=1000)
  r.table=r.table[order(r.table$Fisher,-r.table$Significant),]
  return(list(gdata,r.table))
}

make.GO.table=function(listofGOdata.results,SBIDs,SBGenes,GOnrow=5,GO2Genedb=org.Hs.egGO2ALLEGS){
  allgs=unique(unlist(SBGenes))
  SBSizes=sapply(SBGenes,length)
  SBGenes.c=sapply(SBGenes,function(x)paste(x,collapse=', '))
  GOIDs=c()
  GONms=c()
  GOSizes=c()
  a.An.s=c()
  An.s=c()
  An.s.e=c()
  An.Genes=c()
  GO.ps=c()
  rns=c()
  for(i in 1:length(listofGOdata.results)){
    sbGOdata=listofGOdata.results[[i]][[1]]
    GOt=listofGOdata.results[[i]][[2]]
    t1=min(GOnrow,nrow(GOt))
    rns=c(rns,t1)
    GOt=GOt[1:t1,]
    GOIDs=c(GOIDs,GOt$GO.ID)
    GONms=c(GONms,GOt$Term)
    GOSizes=c(GOSizes,sapply(GOt$GO.ID,function(x){
      count.mappedLkeys(GO2Genedb[x])
    }))
    a.An.s=c(a.An.s,GOt$Annotated)
    An.s=c(An.s,GOt$Significant)
    An.s.e=c(An.s.e,GOt$Expected)
    sbAngs=lapply(GOt$GO.ID,function(x){
      gs=genesInTerm(sbGOdata,x)[[1]]
      return(SBGenes[[i]][SBGenes[[i]] %in% gs])
    })
    An.Genes=c(An.Genes,sapply(sbAngs,function(x)paste(x,collapse=', ')))
    GO.ps=c(GO.ps,GOt$Fisher)
  }
  GO.table=data.frame(rep(SBIDs,rns),rep(SBSizes,rns),rep(SBGenes.c,rns),
                      GOIDs,GONms,GOSizes,a.An.s,An.s,An.s.e,An.Genes,GO.ps)
  return(GO.table)
}

st.title=c("Subnetwork ID","Subnetwork Size","Subnetwork Genes",
           "GO Term ID","GO Term Name","GO Term Size","Network Annotated Size",
           "Subnetwork Annotated Size","Expected Subnetwork Annotated Size",
           "Subnetwork Annotated Genes","p-value")
