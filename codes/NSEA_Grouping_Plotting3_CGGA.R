#Parameters
Folder.Path='TCGA'
Cancer.Code='LGG'
Delete.Data=T #Delete data after loading
par(mar=c(5.1,4.1,4.1,2.1))

#Load data
filepath=paste0(Folder.Path,'/','TCGA_UCSC_',Cancer.Code)
datafilenm=paste0(filepath,'/TCGA_',Cancer.Code,'_Data','.rda')
dataname=paste0('TCGA.',Cancer.Code,'.Data')
load(datafilenm)
save.data=get(dataname)
datafilenm2=paste0(filepath,'/TCGA_NSEA_',Cancer.Code,'_Data','.rda')
dataname2=paste0('TCGA.NSEA.',Cancer.Code,'.Data')
load(datafilenm2)
NSEA.save.data=get(dataname2)

Sbn.Exp=save.data$Sbn.Exp
Sur.l=save.data$Sur.l
Sur.d=save.data$Sur.d
save.data.nms=names(save.data)
if('Subtypes' %in% save.data.nms){
  Subtypes=save.data$Subtypes
}
if('Hist.Type' %in% save.data.nms){
  Hist.Type=save.data$Hist.Type
}
if('Expression.Type' %in% save.data.nms){
  Expression.Type=save.data$Expression.Type
}
if('Methylation.Type' %in% save.data.nms){
  Methylation.Type=save.data$Methylation.Type
}
if('Mutation.Table' %in% save.data.nms){
  Mutation.Table=save.data$Mutation.Table
}
if('Subtype.Table' %in% save.data.nms){
  Subtype.Table=save.data$Subtype.Table
}

Sbn.Network=NSEA.save.data$Sbn.Network
diff.sd=NSEA.save.data$diff.sd
diff.cor.m=NSEA.save.data$diff.cor.m
Sbn.Network3=NSEA.save.data$Sbn.Network3
Sbn.List=NSEA.save.data$Sbn.List
all.sns.gs=NSEA.save.data$all.sns.gs
Sbn.Network21=NSEA.save.data$Sbn.Network21
Sbn.l.cors=NSEA.save.data$Sbn.l.cors
Sbn3.vns=NSEA.save.data$Sbn3.vns
f.Sbn.List=NSEA.save.data$f.Sbn.List
sbn.maxs=NSEA.save.data$sbn.maxs
sbn.maxs2=NSEA.save.data$sbn.maxs2
NSEA.pars=NSEA.save.data$NSEA.pars

Network.file=NSEA.pars$Network.file
Sb.len=NSEA.pars$Sb.len
SD.Threshold=NSEA.pars$SD.Threshold
Gene.Threshold=NSEA.pars$Gene.Threshold
Sbn.Cor.Threshold=NSEA.pars$Sbn.Cor.Threshold
Network.Expansion=NSEA.pars$Network.Expansion
Expansion.Exlusive=NSEA.pars$Expansion.Exlusive

#rm(save.data,NSEA.save.data,NSEA.pars)
#if(Delete.Data){
 # rm(list=c(dataname,dataname2))
#}

#Generate subnetwork matrix
names(f.Sbn.List)=paste0('SN',1:length(f.Sbn.List))
f.Sbn.List1=f.Sbn.List

egv.m=sapply(f.Sbn.List1,function(gs){
  d.m=Sbn.Exp[gs,]
  egv=prcomp(d.m,scale=F,center=F)$rotation[,1]
  egv=remove.outliers(egv,
                  coef1=Outliers.par$coef,
                  a1=Outliers.par$a,
                  b1=Outliers.par$b)
  return(egv)
})

#Consensus Clustering
wardD2 = function(this_dist,k){
  tmp = hclust(this_dist,method='ward.D2')
  assignment = cutree(tmp,k)
  return(assignment)
}

result=ConsensusClusterPlus(t(egv.m),maxK=8,reps=1000,pItem=0.75,pFeature=1,
                            title=paste0(filepath,"/consensus_cluster/"),clusterAlg="pam",distance="euclidean",
                            #innerLinkage="ward.D2", finalLinkage="average",
                            plot="png",seed=1,verbose=T)
egv.memb=result[[5]]$consensusClass
egv.memb.f=factor(egv.memb)
egv.dendro=as.dendrogram(result[[5]]$consensusTree,labels=egv.memb)

#Plot survival curves
Plot.Group=NULL
Surv.Group=factor(egv.memb[names(Sur.l)])

if(is.null(Plot.Group)){
  t1= rep(T,length(Surv.Group))
}else{
  t1= egv.memb %in% Plot.Group
}
Sur.d1=Sur.d[t1]
Sur.l1=Sur.l[t1]
SurvObj=Surv(Sur.d1,event=Sur.l1==1,type='right')
egv.memb.f1=Surv.Group[t1]
sf=survfit(SurvObj ~ egv.memb.f1, conf.type = "log-log")
group.colors=c('#BEBEBE',brewer.pal((max(egv.memb)), 'Set1'))
flvs=as.integer(as.vector(sort(unique(egv.memb.f1))))
survtest.result=survdiff(SurvObj ~ egv.memb.f1,rho=0)
survtest.p.val <- 1 - pchisq(survtest.result$chisq,
                             length(survtest.result$n) - 1)
survtest.p.val=signif(survtest.p.val,2)
survplotnm=paste0(filepath,'/',Cancer.Code,'_Survival_Curves.pdf')
pdf(survplotnm,width=8,height=6)
plot(sf,conf.int=F,col=group.colors[flvs+1],mark.time=T,
     xlab='Time (Months)',ylab='Proportion Survival')
title(paste0(Cancer.Code,' Survival Curves'),line=2.2)
mtext(paste0('p = ',survtest.p.val),side=3,
      line=0.5,font=3,cex=0.8)
lg.text=paste0('Group ',flvs)
lg.text[lg.text=='Group 0']='Ungrouped'
legend(max(Sur.d1),1,lg.text,lty=1,
       col=group.colors[flvs+1],xjust=1,yjust=1)
dev.off()

#Plot heatmap (including patient groups and subtypes)
group.col=group.colors[as.integer(as.vector(egv.memb.f))+1]
colsidecolvector=as.matrix(group.col)
colv.cns='Group'

if('Subtype.Table' %in% save.data.nms){
  mt=Subtype.Table
  mt=as.matrix(mt)
  rownames(mt)=colnames(Sbn.Exp)
  mt[mt=='NA']=NA
  mt.legend=apply(mt,2,make.color)#,qualitative)
  mt.col=sapply(mt.legend,function(x)x$color)
  colsidecolvector=cbind(mt.col,colsidecolvector)
  colv.nnm=colnames(mt.col)
  colv.cns=c(colv.nnm,colv.cns)
}

colnames(colsidecolvector)=colv.cns

s.m=t(egv.m)
t11=quantile(s.m,0.95)
t12=quantile(s.m,0.05)
s.m[s.m<t12]=t12
s.m[s.m>t11]=t11
heatmap.nm=paste0(filepath,'/',Cancer.Code,'_Heatmap.pdf')
pdf(heatmap.nm,width=16,height=8)
heatmap.31(s.m,col=heat.col3,margin=c(5,10),cexRow=0.5,ColSideColors = colsidecolvector,labCol=NA,
           Colv=egv.dendro,lwid=c(1,10),lhei=c(1,2,6),
           key=F)
dev.off()

for(i in 1:length(mt.legend)){
  t1=colnames(mt.col)[i]
  t2=mt.legend[[i]]
  temp.nm=paste0(filepath,'/',Cancer.Code,'_Heatmap_',
                 gsub('\\.','_',t1),
                 '_Legend.pdf')
  pdf(file=temp.nm,width=6,height=4)
  plot.legend(t2,main=t1)
  dev.off()
}


#Network State Plot
fbn3.index=paste0('SN',1:length(f.Sbn.List))
fbn3.index=fbn3.index[fbn3.index %in% names(f.Sbn.List1)]
sn.sta=1
sn.end=length(f.Sbn.List1)
f.Sbn.List3=f.Sbn.List1[fbn3.index[sn.sta:sn.end]]

grps=max(egv.memb)
sub.n=length(f.Sbn.List3)
l.wids=c(1,rep(5,sub.n))
l.hgts=rep(c(1,rep(5,grps)),length(l.wids))
m.c.ns1=names(f.Sbn.List3)
m.c.ns2=paste0('Group',1:grps)
m.c.ns2.col=group.colors[1:grps+1]
l.matrix=matrix(1:((grps+1)*(sub.n+1)),nrow=grps+1)

nsnm=paste0(filepath,'/',Cancer.Code,'_Network_State_',sn.sta,'-',sn.end,'.pdf')
pdf(file=nsnm,width=sum(l.wids),height=sum(l.hgts)/length(l.wids))
#layout(l.matrix,widths=l.wids,heights=l.hgts)
ml1=max(abs(quantile(Sbn.Exp,c(0.05,0.95))))
cmar=par()$mar
par(mar=c(0,0,0,0))
plot.new()
box('plot',lwd=3,col='gray50')
for(j in 1:grps){
  plot.new()
  text(0.5,0.5,m.c.ns2[j],cex=5,srt=90,col=m.c.ns2.col[j])
  box('plot',lwd=3,col='gray50')
}
for(i in 1:sub.n){
  plot.new()
  text(0.5,0.5,m.c.ns1[i],cex=5)
  box('plot',lwd=3,col='gray50')
  gs=f.Sbn.List3[[i]]
  g1=induced.subgraph(Sbn.Network,gs)
  g2=induced.subgraph(Sbn.Network3,gs)
  E(g2)$lty=1
  g1=graph.union(g1,g2)
  E(g1)$lty[is.na(E(g1)$lty)]=2
  g1$layout=layout.auto(g1)
  V(g1)$frame.color='gray75'
  V(g1)$size=25
  V(g1)$label.font=2
  V(g1)$label.cex=2
  V(g1)$label.color='black'
  for(j in 1:grps){
    v.cols=get.color2(V(g1)$name,egv.memb==j,Sbn.Exp,ml1)
    plot(g1,vertex.color=v.cols)
    # plot(g1,vertex.size=V(g1)$size*2.5)
    box('plot',lwd=3,col='gray50')
  }
}
par(mar=cmar)
dev.off()

#Gene set enrichment analysis
Selected.Subnetwork=NULL
Gene.Pool=rownames(Sbn.Exp)
Geneset.Size.min=1000
GSEA.Algorithm='classic'
GSEA.Statistic='fisher'
Mapping.Rule='org.Hs.eg.db'
GO2Gene.db='org.Hs.egGO2ALLEGS'
GO.N=10

if(file.exists('sbngsea_parallel.log')){
  file.remove('sbngsea_parallel.log')
}
if(exists('CLS')){
  try(stopCluster(CLS))
}
closeAllConnections()
CLS=makeCluster(10,outfile='sbngsea_parallel.log')
registerDoParallel(CLS)

if(is.null(Selected.Subnetwork)){
  Select.Sbn=1:length(f.Sbn.List3)
}else{
  Select.Sbn=Selected.Subnetwork
}
gs=NULL
allgs=Gene.Pool
ns=Geneset.Size.min
alg=GSEA.Algorithm
teststa=GSEA.Statistic
mapping.base=Mapping.Rule

ags=rep(0,length(allgs))
names(ags)=allgs
gdata0=data=new("topGOdata",description='Subnetwork Genes',ontology="BP",allGenes=ags,
                geneSel=function(x)x==1,nodeSize=ns,annot=annFUN.org,mapping=mapping.base,ID="symbol")

clusterExport(CLS,c('gdata0','allgs','GOtest','alg','teststa'))
Current.sub.GSEA=parLapply(CLS,f.Sbn.List3[Select.Sbn],function(gs){
  library(topGO)
  gdata=gdata0
  ags=rep(0,length(allgs))
  names(ags)=allgs
  ags[allgs %in% gs]=1
  gdata=updateGenes(gdata,ags,function(x)x==1)
  result=GOtest(gdata,alg,teststa)
  r.table=GenTable(gdata,Fisher=result,topNodes=20,numChar=1000)
  r.table=r.table[order(r.table$Fisher,-r.table$Significant),]
  return(list(gdata,r.table))
})

t16=lapply(Current.sub.GSEA,function(x)x[[2]])
t17=sapply(t16,function(x)x$Term[[1]])
t17.1=paste0(as.character(m.c.ns1[Select.Sbn]),': ',t17)
t18=sapply(t16,function(x)x$Fisher[[1]])
t19=sapply(f.Sbn.List3[Select.Sbn],length)
t110=sapply(t16,function(x)x$Significant[[1]])
Current.sub.stat=data.frame(ID=Select.Sbn,SBSize=t19,BF=t17.1,SBSig=t110,Pv=t18)
Current.sub.stat$Unsig=Current.sub.stat$SBSize-Current.sub.stat$SBSig

GO2Gene.db1=get(GO2Gene.db)
Current.sub.stat.t=make.GO.table(Current.sub.GSEA,m.c.ns1[Select.Sbn],
                                 f.Sbn.List3[Select.Sbn],GOnrow=GO.N,GO2Genedb=GO2Gene.db1)
Current.sub.stat.t=rbind(st.title,Current.sub.stat.t)
rownames(Current.sub.stat.t)=NULL
colnames(Current.sub.stat.t)=paste0('V',1:ncol(Current.sub.stat.t))

gsea.t.nm=paste0(filepath,'/',Cancer.Code,'_GSEA_Table.tab')
export.2(Current.sub.stat.t,gsea.t.nm)

bar.col='black'
bar.mn.max=90

bar.m=t(as.matrix(Current.sub.stat[,c('SBSig','Unsig')]))[,nrow(Current.sub.stat):1]
bar.mn=substr(rev(t17.1),1,bar.mn.max)
bar.l=rev(t18)
bar.xm=ceiling(max(Current.sub.stat$SBSize)/2)*2

gsea.bar.nm=paste0(filepath,'/',Cancer.Code,'_GSEA.pdf')
pdf(file=gsea.bar.nm,height=64,width=16)
cmar=par()$mar
par(mar=cmar+c(0,35,0,15))
par(xpd=T)
bp=barplot(bar.m,density=0,border=0,
           las=2,horiz=T,cex.names=1,names.arg=bar.mn,
           xaxt="n",
           xlab='Number of Genes',main='RNA Subnetwork Function',cex.main=2)
axis(1, at = seq(0, bar.xm, by = 2))
legend('right', inset=c(-0.33,0), fill=T,legend=c('Annotated','Unannotated'),cex=1.3, bty='n',col=1, angle=45, density=10)
legend('right', inset=c(-0.33,0), fill=T,legend=c('Annotated','Unannotated'),cex=1.3, bty='n',col=1, angle=135, density=c(10,0))
text(x=colSums(bar.m),y=bp,labels=bar.l,pos=4)
barplot(bar.m,
        las=2,horiz=T,cex.names=1,col=bar.col,density=10,axes=F,angle=45,names.arg=rep(NA,ncol(bar.m)),
        xlim=c(0,bar.xm),add=T)
barplot(bar.m,
        las=2,horiz=T,cex.names=1,col=bar.col,density=c(10,0),angle=135,border=F,axes=F,
        names.arg=rep(NA,ncol(bar.m)),xlim=c(0,bar.xm),add=T)
par(mar=cmar)
par(xpd=F)
dev.off()
