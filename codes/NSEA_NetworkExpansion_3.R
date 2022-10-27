Sbn.Network22=Sbn.Network21

gnr=function(g){
  d1=graph.bfs(g, root=1, 
               order=F, rank=F, father=F, pred=F,
               succ=F, dist=T)$dist
  t1=d1 %% 2==0
  x=sum(t1)/sum(!t1)
  if(x<1) x=1/x
  return(x)
}

if(Expansion.Exlusive){
  f.Sbn.List=list()
  r.l=Sbn.List
  s.s=Sbn.l.cors
  count=0
  sbn.maxs=c()
  sbn.maxs2=c()
  repeat{
    count=count+1
    t1=which.max(s.s)
    sbns.max=max(s.s)
    sbn.maxs2=c(sbn.maxs2,sbns.max)
    max.gs=r.l[,t1]
    
    #--------------------------------
    c.c.gs=max.gs
    c.neighbors=unlist(lapply(c.c.gs,function(g)neighbors(Sbn.Network22,g)$name))
    c.neighbors=c.neighbors[!c.neighbors %in% c.c.gs]
    gs=c.c.gs
    gs=Sbn3.vns[Sbn3.vns %in% gs]
    m=Sbn.Adj.m[gs,gs]
    m[lower.tri(m,F)]=0
    c.il=which(m==1,arr.ind=T)
    c.gl=data.table(apply(c.il,2,function(x)gs[x]))
    eg.is=egl21[c.gl]$eg.ind
    corm1=diff.cor.m[eg.is,eg.is]
    corm1=corm1[lower.tri(corm1)]
    p.diff=mean(abs(corm1))
    p.degree=gnr(induced_subgraph(Sbn.Network22,c.c.gs))
    # p.degree=edge_density(induced_subgraph(Sbn.Network22,c.c.gs))
    # p.dxd=p.diff+p.degree
    repeat{
      if(length(c.neighbors)==0) break
      results=sapply(c.neighbors,function(gene){
        gs=c(c.c.gs,gene)
        gs=Sbn3.vns[Sbn3.vns %in% gs]
        m=Sbn.Adj.m[gs,gs]
        m[lower.tri(m,F)]=0
        c.il=which(m==1,arr.ind=T)
        c.gl=data.table(apply(c.il,2,function(x)gs[x]))
        eg.is=egl21[c.gl]$eg.ind
        corm1=diff.cor.m[eg.is,eg.is]
        corm1=corm1[lower.tri(corm1)]
        diff=mean(abs(corm1))
        g1=induced_subgraph(Sbn.Network22,gs)
        d1=gnr(g1)
        return(c(diff,d1))
      })
      results=cbind(results)
      diffs=results[1,]
      degrees=results[2,]
      # dxd=diffs+degrees
      t1=(diffs>p.diff | diffs>0.8) & degrees<=2
      if(any(t1)){
        # dxd1=dxd[t1]
        #     if(max(dxd1)<p.dxd){
        #       break
        #     }
        diffs1=diffs[t1]
        degrees1=degrees[t1]
        c.neighbors1=c.neighbors[t1]
        t2=which.max(diffs1)
        # t2=which.min(degrees1)
        t.g=c.neighbors1[t2]
        # print(t.g)
        print(diffs1[t2])
        p.diff=diffs1[t2]
        p.degree=degrees1[t2]
        # p.dxd=max(dxd1)
        # print(p.diff)
        # print(p.degree)
        c.c.gs=c(c.c.gs,t.g)
        n.gs=neighbors(Sbn.Network22,t.g)$name
        c.neighbors=unique(c.neighbors,n.gs)
        c.neighbors=c.neighbors[!c.neighbors %in% c.c.gs]
      }else{
        break
      }
    }
    if(length(c.c.gs)>=5){
      max.gs=c.c.gs
      plot(induced_subgraph(Sbn.Network22,max.gs))
      Sbn.Network22=Sbn.Network22-max.gs
      names(max.gs)=NULL
      sbns.max=p.diff
      print(max.gs)
      sbn.maxs=c(sbn.maxs,sbns.max)
      f.Sbn.List=c(f.Sbn.List,list(max.gs))
      r.l.l=apply(r.l,1,function(x)x %in% max.gs)
      if(is.null(dim(r.l.l))) r.l.l=t(r.l.l)
      # t2=rowSums(r.l.l)<=floor(length(max.gs)*0.4)
      t2= rowSums(r.l.l)==0
      if(!any(t2)) break
      r.l=r.l[,t2,drop=F]
      s.s=s.s[t2]
    }else{
      print(max.gs)
      if(NCOL(r.l)==1) break
      t1=which.max(s.s)
      r.l=r.l[,-t1,drop=F]
      s.s=s.s[-t1]
    }
  }
  # names(f.Sbn.List)=paste0('V',1:length(f.Sbn.List))
}else{
  # f.Sbn.Lists=list()
  f.Sbn.List=list()
  #   c.cth.i=1
  #   c.cth=cths[c.cth.i]
  r.l=Sbn.List
  s.s=Sbn.l.cors
  count=0
  sbn.maxs=c()
  repeat{
    count=count+1
    t1=which.max(s.s)
    sbns.max=max(s.s)
    #     while(sbns.max<c.cth){
    #       f.Sbn.Lists=c(f.Sbn.Lists,list(f.Sbn.List))
    #       c.cth.i=c.cth.i+1
    #       if(c.cth.i>length(cths)) break
    #       c.cth=cths[c.cth.i]
    #     }
    #     if(c.cth.i>length(cths)) break
    max.gs=r.l[,t1]
    names(max.gs)=NULL
    print(max.gs)
    sbn.maxs=c(sbn.maxs,sbns.max)
    f.Sbn.List=c(f.Sbn.List,list(max.gs))
    r.l.l=apply(r.l,1,function(x)x %in% max.gs)
    if(is.null(dim(r.l.l))) r.l.l=t(r.l.l)
    t2= rowSums(r.l.l)==0
    if(!any(t2)) break
    r.l=r.l[,t2,drop=F]
    s.s=s.s[t2]
  }
  sbn.maxs2=sbn.maxs
  
  clusterExport(CLS,c('Sbn.Network22','Sbn.Adj.m','egl21',
                      'diff.cor.m','Sbn3.vns','gnr'))
  results=parLapply(CLS,f.Sbn.List,function(c.c.gs){
    library(igraph)
    library(data.table)
    c.neighbors=unlist(lapply(c.c.gs,function(g)neighbors(Sbn.Network22,g)$name))
    c.neighbors=c.neighbors[!c.neighbors %in% c.c.gs]
    gs=c.c.gs
    gs=Sbn3.vns[Sbn3.vns %in% gs]
    m=Sbn.Adj.m[gs,gs]
    m[lower.tri(m,F)]=0
    c.il=which(m==1,arr.ind=T)
    c.gl=data.table(apply(c.il,2,function(x)gs[x]))
    eg.is=egl21[c.gl]$eg.ind
    corm1=diff.cor.m[eg.is,eg.is]
    corm1=corm1[lower.tri(corm1)]
    p.diff=mean(abs(corm1))
    p.degree=gnr(induced_subgraph(Sbn.Network22,c.c.gs))
    repeat{
      if(length(c.neighbors)==0) break
      results=sapply(c.neighbors,function(gene){
        gs=c(c.c.gs,gene)
        gs=Sbn3.vns[Sbn3.vns %in% gs]
        m=Sbn.Adj.m[gs,gs]
        m[lower.tri(m,F)]=0
        c.il=which(m==1,arr.ind=T)
        c.gl=data.table(apply(c.il,2,function(x)gs[x]))
        eg.is=egl21[c.gl]$eg.ind
        corm1=diff.cor.m[eg.is,eg.is]
        corm1=corm1[lower.tri(corm1)]
        diff=mean(abs(corm1))
        g1=induced_subgraph(Sbn.Network22,gs)
        d1=gnr(g1)
        return(c(diff,d1))
      })
      results=cbind(results)
      diffs=results[1,]
      degrees=results[2,]
      # dxd=diffs+degrees
      t1=diffs>p.diff & degrees<=2
      if(any(t1)){
        # dxd1=dxd[t1]
        #     if(max(dxd1)<p.dxd){
        #       break
        #     }
        diffs1=diffs[t1]
        degrees1=degrees[t1]
        c.neighbors1=c.neighbors[t1]
        t2=which.max(diffs1)
        # t2=which.min(degrees1)
        t.g=c.neighbors1[t2]
        # print(t.g)
        p.diff=diffs1[t2]
        p.degree=degrees1[t2]
        # p.dxd=max(dxd1)
        # print(p.diff)
        # print(p.degree)
        c.c.gs=c(c.c.gs,t.g)
        n.gs=neighbors(Sbn.Network22,t.g)$name
        c.neighbors=unique(c.neighbors,n.gs)
        c.neighbors=c.neighbors[!c.neighbors %in% c.c.gs]
      }else{
        break
      }
    }
#     repeat{
#       if(length(c.neighbors)==0) break
#       results=sapply(c.neighbors,function(gene){
#         gs=c(c.c.gs,gene)
#         gs=Sbn3.vns[Sbn3.vns %in% gs]
#         m=Sbn.Adj.m[gs,gs]
#         m[lower.tri(m,F)]=0
#         c.il=which(m==1,arr.ind=T)
#         c.gl=data.table(apply(c.il,2,function(x)gs[x]))
#         eg.is=egl21[c.gl]$eg.ind
##         diff=mean(abs(diff.cor.m[eg.is,eg.is]))
#         return(diff)
#       })
#       if(max(results)>p.diff){
#         t.g=c.neighbors[which.max(results)]
#         # print(t.g)
#         p.diff=max(results)
#         # print(p.diff)
#         c.c.gs=c(c.c.gs,t.g)
#         n.gs=neighbors(Sbn.Network22,t.g)$name
#         c.neighbors=unique(c.neighbors,n.gs)
#         c.neighbors=c.neighbors[!c.neighbors %in% c.c.gs]
#       }else{
#         break
#       }
#     }
    max.gs=c.c.gs
    names(max.gs)=NULL
    sbns.max=p.diff
    # print(max.gs)
    return(list(c.c.gs=c.c.gs,sbns.max=sbns.max))
  })
  f.Sbn.List=lapply(results,function(x)x[[1]])
  sbn.maxs=sapply(results,function(x)x[[2]])
}

#save('f.Sbn.List','sbn.maxs','sbn.maxs2',file='degree_limit_ne')


