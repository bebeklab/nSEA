hist.type=pat.m[,'Histology']
names(hist.type)=make.names(pat.m[,'Patient'])

# fp1=paste0(filepath,'/','LGG_subtypes.txt')
# Subtypes.LGG=import.1(fp1,T)
# Sbts=Subtypes.LGG[,2]
# names(Sbts)=make.names(Subtypes.LGG[,1])
# Sbts[Sbts=='NA']=NA
# mutation.table=Subtypes.LGG[,4:ncol(Subtypes.LGG)]
# rownames(mutation.table)=make.names(Subtypes.LGG[,1])

fp1=paste0(filepath,'/','tcgaclusters.txt')
Subtypes.LGG2=import.1(fp1,T)
# Sbts2=Subtypes.LGG2[,2]
# names(Sbts2)=make.names(Subtypes.LGG2[,1])
# Sbts2[Sbts2=='NA']=NA
subtype.table=Subtypes.LGG2[,2:ncol(Subtypes.LGG2)]
rownames(subtype.table)=make.names(Subtypes.LGG2[,1])

hist.type=hist.type[colnames(save.data$Sbn.Exp)]
subtype.table=subtype.table[colnames(save.data$Sbn.Exp),]
# Sbts=Sbts[colnames(save.data$Sbn.Exp)]
# mutation.table=mutation.table[colnames(save.data$Sbn.Exp),]
# colnames(mutation.table)[1]='TERT.promoter'
# save.data$Subtypes=Sbts
save.data$Hist.Type=hist.type
save.data$Subtype.Table=subtype.table
# save.data$Mutation.Table=mutation.table


# Sbts2=Sbts2[colnames(save.data$Sbn.Exp)]

# save.data$Subtypes2=Sbts2
# save.data$Hist.Type=hist.type

save.name=paste0('TCGA.',Cancer.Code,'.Data')
save.name1=paste0(filepath,'/TCGA_',Cancer.Code,'_Data','.rda')
assign(save.name,save.data,envir=.GlobalEnv)
save(list=save.name,file=save.name1)
#rm(save.data)
save.data.par=list(Transformation=Transformation,
                   Outliers.Rm=Outliers.Rm,
                   Outliers.par=Outliers.par)
save.name21=paste0('TCGA.',Cancer.Code,'.Data.par')
save.name2=paste0(filepath,'/TCGA_',Cancer.Code,'_Data_par','.rda')
assign(save.name21,save.data.par,envir=.GlobalEnv)
save(list=save.name21,file=save.name2)
