# This script counts ribosomal genes in the DE and DS genes found for each dataset, and plots figure 5 using the results.
# Its input are the HBA-DEALS output files from the Snakemake pipeline.

library(biomaRt)
library(ggpubr)

source('get_fdr_prob.R')

ribo.mat=matrix(nrow=0,ncol=3)

colnames(ribo.mat)=c('dataset','DAS','DGE')

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host="www.ensembl.org")


for (dataset in c("mason_latest",
                  "SRP040070_3","SRP040070_7","SRP040070_9",
                  "SRP227272_38",'SRP279203_72','SRP294125_74',"SRP078309_53","SRP178454_50","SRP186406_51","SRP216763_55","SRP222569_54","SRP251704_52",
                  "SRP273785_56",'SRP284977_76','SRP284977_77'))
{
  
  res=read.table(paste0('HBA-DEALS-covid19_output/',dataset,'.txt'),header=TRUE)
  
  gene.names=getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                   filters    = "ensembl_gene_id",
                   values     = res$Gene, 
                   mart       = mart)
  
  gene.labels=unlist(lapply(res$Gene,function(g)gene.names$hgnc_symbol[which(gene.names$ensembl_gene_id==g)][1]))
  
  if (dataset=='mason_latest')
    
    gene.labels=res$Gene
  
  thresh=get.fdr.prob(res)
  
  de.genes=gene.labels[res$Isoform=='Expression' & res$P<=thresh[[1]]]
  
  ds.genes=unique(gene.labels[res$Isoform!='Expression' & res$P<=thresh[[2]]])
  
  ribo.mat=rbind(ribo.mat,c(dataset,sum(grepl('RPS',ds.genes) | grepl('RPL',ds.genes))/length(ds.genes),sum(grepl('RPS',de.genes) | grepl('RPL',de.genes))/length(de.genes)))
  
}

write.table(ribo.mat,'RPS_genes_fraction.txt',sep='\t',quote = F,row.names = F,col.names = T)


df=data.frame(fraction=as.numeric(as.character(ribo.mat[1:16,2])),dataset=ribo.mat[1:16,1],type=ifelse(1:16 %in% c(1:7,15:16),'Coronavirus','Other'))


sarscov2=c('mason_latest','SRP279203_72','SRP294125_74','SRP284977_76','SRP284977_77')

df$category <- as.factor(df$dataset %in% sarscov2)

g <- ggplot(df, aes(x=type, y=fraction)) +
  geom_boxplot(alpha = 0.05,colour="black",fill="blue") +
  geom_point(aes(shape=category,color=category,size=2),position = position_jitterdodge()) +
  scale_color_manual(values=c("#E64B36", "#3C5488")) +
  theme_linedraw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("")+ theme(axis.text = element_text(size = 20))+ theme(axis.title.y  = element_text(size = 20)) 

g

wilcox.test(df$fraction[c(1:7,15:16)],df$frac[8:14],paired=FALSE)

df=data.frame(fraction=as.numeric(as.character(ribo.mat[1:16,3])),dataset=ribo.mat[1:16,1],type=ifelse(1:16 %in% c(1:7,15:16),'Coronavirus','Other'))

df$category <- as.factor(df$dataset %in% sarscov2)

g <- ggplot(df, aes(x=type, y=fraction)) +
  geom_boxplot(alpha = 0.05,colour="black",fill="blue") +
  geom_point(aes(shape=category,color=category,size=2),position = position_jitterdodge()) +
  scale_color_manual(values=c("#E64B36", "#3C5488")) +
  theme_linedraw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("")+ theme(axis.text = element_text(size = 20))+ theme(axis.title.y  = element_text(size = 20)) 


g

wilcox.test(df$frac[c(1:7,15:16)],df$frac[8:14],paired=FALSE)


