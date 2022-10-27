#Parameters
Folder.Path='TCGA'
Cancer.Code='LGG'

Network.file='H.GGI.string' #PPI network file * .7
SD.Threshold=0.95 #Edge pre-filtering threshold

#Fixed parameters
Sb.len=4 #Subnetwork size
Network.Expansion=T #Whether to expand the subnetwork features
Expansion.Exlusive=T #Whether the expansion is exclusive

#Load the data file
filepath=paste0(Folder.Path,'/','TCGA_UCSC_',Cancer.Code)
datafilenm=paste0(filepath,'/TCGA_',Cancer.Code,'_Data','.rda')
dataname=paste0('TCGA.',Cancer.Code,'.Data')
load(datafilenm)
save.data=get(dataname)

Sbn.Exp=save.data$Sbn.Exp # expression matrix
Sur.l=save.data$Sur.l
Sur.d=save.data$Sur.d

#Set up parallel processing
if(file.exists('sbncal_parallel.log')){
  file.remove('sbncal_parallel.log')
}
if(exists('CLS')){
  try(stopCluster(CLS))
}
closeAllConnections()

CLS=makeCluster(10,outfile='sbncal_parallel.log')
registerDoParallel(CLS)

#Load PPI network and filter out genes not included in the expression matrix
network.nm=load(Network.file)
H.GGI.string=get(network.nm)
Sbn.Network=induced.subgraph(H.GGI.string,
                             intersect(rownames(Sbn.Exp),V(H.GGI.string)$name))
Sbn.Network=remove.graph.islands(Sbn.Network,min.v=Sb.len)

#Edge filtration
t1=apply(Sbn.Exp,1,function(x)length(unique(x)))
egl=get.edgelist(Sbn.Network)
clusterExport(CLS,c('egl','Sbn.Exp'))
diff.sd=parSapply(CLS,1:nrow(egl),function(i){
  x=egl[i,]
  t1=Sbn.Exp[x[1],]
  t2=Sbn.Exp[x[2],]
  sd1=sd(t1-t2)
  return(sd1)
})

t1=diff.sd>quantile(diff.sd,SD.Threshold)
egl2=egl[t1,]
diff.sd2=diff.sd[t1]
Sbn.Network2=Sbn.Network-E(Sbn.Network)[!t1]

#Edge correlation matrix
clusterExport(CLS,c('egl2','diff.sd2'))
diff.m2=parSapply(CLS,1:nrow(egl2),function(i){
  x=egl2[i,]
  t1=Sbn.Exp[x[1],]
  t2=Sbn.Exp[x[2],]
  return(t1-t2)
})
diff.cor.m=cor(diff.m2)
diag(diff.cor.m)=0

#Enumeration
Sbn.Network3=remove.graph.islands(Sbn.Network2,Sb.len)
Sbn.Network21=induced.subgraph(Sbn.Network2,V(Sbn.Network3)$name)

Decomposition4.3=function(d.graph){
  n3.list=pblapply(V(d.graph)$name,function(x){
    gs=names(neighbors(d.graph,x))
    if(length(gs)>=2){
      return(rbind(x,as.matrix(combn(gs,2))))
    }else{
      return(character())
    }
  })
  n3.list1=do.call(cbind,n3.list)
  n3.list1=apply(n3.list1,2,sort)
  n3.list1=t(unique(t(n3.list1)))
  return(n3.list1)
}

#Generate 3-node subnetworks
n3.sbn.list=Decomposition4.3(Sbn.Network3)

#Generate 4-node subnetworks
clusterExport(CLS,c('n3.sbn.list','Sbn.Network3'))
n4.sbn.list.raw=parLapply(CLS,1:ncol(n3.sbn.list),function(i){
  library(igraph)
  library(data.table)
  if(i %% 1000==0) print(i)
  x=n3.sbn.list[,i]
  gs=lapply(adjacent_vertices(Sbn.Network3,x),names)
  gs=unique(unlist(gs))
  gs=gs[!gs %in% x]
  n.gs=gs
  if(length(n.gs)<1) return(c())
  if(length(n.gs)==1){
    t1=sort(c(x,n.gs))
    names(t1)=paste0('V',1:NROW(t1))
    t1=data.table(t(t1))
    return(t1)
  }
  if(length(n.gs)>1){
    x.gs=rep(x,length(n.gs))
    dim(x.gs)=c(3,length(n.gs))
    cb.gs=rbind(x.gs,n.gs)
    cb.gs=t(apply(cb.gs,2,sort))
    t1=cb.gs
    colnames(t1)=paste0('V',1:ncol(t1))
    t1=data.table(t1)
    return(t1)
  }
})

#Remove duplicated subnetworks
n4.sbn.list=rbindlist(n4.sbn.list.raw)
setkey(n4.sbn.list)
n4.sbn.list=unique(n4.sbn.list)

Sbn.List=as.matrix(n4.sbn.list)
Sbn.List=t(Sbn.List)
all.sns.gs=unique(as.vector(Sbn.List))

#Calculate inner-pattern consistency (subnetwork score)
egl21=egl2
egl21=data.table(egl21)
egl21[,eg.ind:=1:nrow(egl21)]
setkey(egl21,'V1','V2')
Sbn.Adj.m=as_adjacency_matrix(Sbn.Network21,attr=NULL,names=T,edges=F,sparse=F)

Sbn3.vns=V(Sbn.Network21)$name
clusterExport(CLS,c('Sbn3.vns','egl21','Sbn.Adj.m','data.table','diff.cor.m'))
Sbn.l.cors=parApply(CLS,Sbn.List,2,function(gs){
  gs=Sbn3.vns[Sbn3.vns %in% gs]
  m=Sbn.Adj.m[gs,gs]
  ds=colSums(m)
  t1=sum(ds==1)
  if(t1>=(nrow(m)-1)){
    return(0)
  }
  m[lower.tri(m,F)]=0
  c.il=which(m==1,arr.ind=T)
  c.gl=data.table(apply(c.il,2,function(x)gs[x]))
  eg.is=egl21[c.gl]$eg.ind
  corm1=diff.cor.m[eg.is,eg.is]
  corm1=corm1[lower.tri(corm1)]
  return(mean(abs(corm1)))
})

#Network expansion
source('NSEA_NetworkExpansion_3.R')

#Save results to files
NSEA.pars=list(Network.file=Network.file,
               Sb.len=Sb.len,SD.Threshold=SD.Threshold,
               #Gene.Threshold=Gene.Threshold,
               Network.Expansion=Network.Expansion,
               Expansion.Exlusive=Expansion.Exlusive)

NSEA.save.data=list(Sbn.Network=Sbn.Network,diff.sd=diff.sd,
                    diff.cor.m=diff.cor.m,Sbn.Network3=Sbn.Network3,
                    Sbn.List=Sbn.List,
                    all.sns.gs=all.sns.gs,
                    Sbn.Network21=Sbn.Network21,Sbn.l.cors=Sbn.l.cors,
                    Sbn3.vns=Sbn3.vns,f.Sbn.List=f.Sbn.List,
                    sbn.maxs2=sbn.maxs2,sbn.maxs=sbn.maxs,
                    NSEA.pars=NSEA.pars)

save.name=paste0('TCGA.NSEA.',Cancer.Code,'.Data')
save.name1=paste0(filepath,'/TCGA_NSEA_',Cancer.Code,'_Data','.rda')
assign(save.name,NSEA.save.data,envir=.GlobalEnv)
save(list=save.name,file=save.name1)
#rm(NSEA.save.data,NSEA.pars)
