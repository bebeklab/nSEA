#--------------Parameters-------------------------------------------
Folder.Path='TCGA' #Directory
Cancer.Code='LGG' #Cancer type
Transformation='none' #Data transformation process. Options: log2, log2.1, none. log2.1 is log(x+1,2)
Outliers.Rm=T #Whether to remove outliers by adjboxStats function
Outliers.par=list(coef=1.5,a=-4,b=3) #Outlier removing parameters

#---------------Data Processing-----------------

#Extract expression matrix and do data transformation
filepath=paste0(Folder.Path,'/','TCGA_UCSC_',Cancer.Code)
fp=paste0(filepath,'/','genomicMatrix')
exp.d=import.1(fp,T)
rownames(exp.d)=exp.d[,1]
exp.d=exp.d[,2:ncol(exp.d)]
barcodes=colnames(exp.d)
sample.type=substr(barcodes,14,15)
exp.d=exp.d[,sample.type=='01']
barcodes=barcodes[sample.type=='01']
bar.p=substr(barcodes,1,12)
patients=unique(bar.p)
if(Transformation=='log2'){
  exp.d=log(exp.d,2)
}else if(Transformation=='log2.1'){
  exp.d=log(exp.d+1,2)
}
exp.m2=sapply(patients,function(x)rowMeans(exp.d[,bar.p==x,drop=F]))
colnames(exp.m2)=make.names(patients)
exp.m3=exp.m2
exp.m4=exp.m3
rm(exp.m2,exp.d)

#Input clinical data
c.f=paste0(filepath,'/','clinicaldata.txt')
pat.m=import.1(c.f,T)
c.bar=pat.m[,'Patient']
c.bar=make.names(c.bar)
os.t=pat.m[,6]
os.l=pat.m[,7]
c.bar0=c.bar
os.t[os.t=='NA']=NA
os.t=as.numeric(os.t)
os.l[os.l=='NA']=NA
os.l=as.numeric(os.l)
t1=is.na(os.t) | is.na(os.l)
Sur.Data=cbind(os.t,os.l)
colnames(Sur.Data)=c("Survival","Censor")
rownames(Sur.Data)=c.bar
Sur.Data=Sur.Data[!t1,]
Sur.Data=Sur.Data[Sur.Data[,1]!=0,]
t2=colnames(exp.m4)
t2=t2[t2 %in% rownames(Sur.Data)]
clin.d=Sur.Data[t2,]
Sur.d=clin.d[,1]
Sur.l=clin.d[,2]

#Outlier cutoff values
if(Outliers.Rm){
  Sbn.Exp=t(pbapply(exp.m4,1,function(x){
    x1=remove.outliers(x,
                       coef1=Outliers.par$coef,
                       a1=Outliers.par$a,
                       b1=Outliers.par$b)
    return(getZ2(x1))
  }))
}else{
  Sbn.Exp=t(pbapply(exp.m4,1,getZ2))
}
colnames(Sbn.Exp)=colnames(exp.m4)
rm(exp.m3)

#Input additional (clinical data) and save all processed data to a file
save.data=list(Sbn.Exp=Sbn.Exp,Sur.d=Sur.d,Sur.l=Sur.l)
rm(Sbn.Exp)
fn1=paste0('Data_Additional_Processing_',Cancer.Code,'.R')
if(fn1 %in% list.files(filepath)){
  fp1=paste0(filepath,'/',fn1)
  source(fp1)
}
