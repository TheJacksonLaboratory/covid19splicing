library("ggsci")
library("ggplot2")
library("gridExtra")

source('get_fdr_prob.R')

res=read.table(paste0('HBA-DEALS-covid19_output/mason_latest.txt'),header=TRUE)

gene.labels=res$Gene

thresh=get.fdr.prob(res)

ds.genes=unique(gene.labels[res$Isoform!='Expression' & res$P<=thresh[[2]]])

mason.counts=read.table('mason_rsem_matrix.txt',header=TRUE,sep='\t')

counts.to.proportions=function(isoform.counts){
  
  return(do.call(rbind,lapply(split(isoform.counts,isoform.counts[,1]),function(m){
    v=colSums(m[,3:ncol(m)]);m[,c(F,F,v==0)]=1/nrow(m);m[,c(F,F,v!=0)]=t(t(m[,c(F,F,v!=0)])/colSums(m[,c(F,F,v!=0)]));return(m)})))
  
}

mason.counts[,1]=as.character(mason.counts[,1])

mason.counts=counts.to.proportions(mason.counts)

meta.data=readxl::read_xlsx('compare_old_new_metadata.xlsx')

meta.data=meta.data[meta.data$SampleID_NEWMETADATA_23Jun %in% colnames(mason.counts),]

meta.data=meta.data[match(colnames(mason.counts)[-c(1,2)],meta.data$SampleID_NEWMETADATA_23Jun),]

mason.counts=mason.counts[,c(1,2,which(meta.data$Type_NEWMETADATA_23Jun!='OtherViralInfection')+2)]

meta.data=meta.data[meta.data$Type_NEWMETADATA_23Jun!='OtherViralInfection',]

sum(colnames(mason.counts)[-c(1,2)]==meta.data$SampleID_NEWMETADATA_23Jun)==ncol(mason.counts)-2

iso.virload.cor=apply(mason.counts[,-c(1,2)],1,function(x)cor.test(x=x,y=meta.data$SARS_CoV2_fraction_NEWMETADATA_23Jun,method='kendall')$p.value)

df=data.frame(minus.log10.p=-log10(iso.virload.cor),is.das=mason.counts$gene_id %in% ds.genes)

library(ggpubr)

ggboxplot(data=df,x='is.das',y='minus.log10.p',palette = 'npg',fill='is.das')+geom_hline(yintercept=-log10(0.05), linetype='dotted', col = 'red')+
  annotate("text", x = 1.5, y = -log10(0.05), label = "-log10(0.05)", vjust = -0.5,cex=2)
  
  