get.fdr.prob=function(res,fdr=0.05,resolution=0.01)
{
  
  expected.h0.de=unlist(lapply(seq(resolution,0.25,resolution),function(p)sum(res$P[res$Isoform=='Expression'& res$P<=p])/sum(res$Isoform=='Expression'& res$P<=p)))
  
  exp.thresh=seq(resolution,0.25,resolution)[max(which(expected.h0.de<=fdr))]
  
  expected.h0.ds=unlist(lapply(seq(resolution,0.25,resolution),function(p)
    
  sum((res$P)[res$Isoform!='Expression'& res$P<=p])/sum(res$Isoform!='Expression' & res$P<=p)))
  
  spl.thresh=seq(resolution,0.25,resolution)[max(which(expected.h0.ds<=fdr))]
  
  return(list(expression=exp.thresh,splicing=spl.thresh))
  
}