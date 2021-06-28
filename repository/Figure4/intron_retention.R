# This script finds the proportion of retained intron isoforms among the differentially spliced and differentially expressed 
# isoforms of each dataset and plots the result as figure 4.  Its inputs are the HBA-DEALS output files from the Snakemake pipeline.


library(biomaRt)
library(ggpubr)

source('get_fdr_prob.R')

results.mat=matrix(ncol=2,nrow=0)

colnames(results.mat)=c('dataset','frac.ri')

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host="www.ensembl.org")

for (dataset in c("mason_latest",
  "SRP040070_3","SRP040070_7","SRP040070_9",
  "SRP227272_38",'SRP279203_72','SRP294125_74',"SRP078309_53","SRP178454_50","SRP186406_51","SRP216763_55","SRP222569_54","SRP251704_52",
  "SRP273785_56",'SRP284977_76','SRP284977_77'))
{
  
  res=read.table(paste0('HBA-DEALS-covid19_output/',dataset,'.txt'),header=TRUE)
  

  isoform.types=getBM(attributes = c("ensembl_transcript_id","transcript_biotype"),
                  filters    = "ensembl_transcript_id",
                  values     = res$Isoform[res$Isoform!='Expression'] , 
                  mart       = mart)
  
  res=res[res$Isoform=='Expression' | res$Isoform %in% isoform.types$ensembl_transcript_id,]
  
  isoform.types=isoform.types[match(res$Isoform[res$Isoform!='Expression'],isoform.types$ensembl_transcript_id),]
  
  res=cbind(res,rep('Gene',nrow(res)))
  
  colnames(res)[ncol(res)]='Type'
  
  res$Type=as.character(res$Type)
  
  res$Type[res$Isoform!='Expression']=as.character(isoform.types$transcript_biotype)
  
  thresh=get.fdr.prob(res)
  
  ds.isoforms=res$Isoform[res$Isoform!='Expression' & res$P<=thresh[[2]]]
  
  types.tab=table(res$Type[(res$Isoform %in% ds.isoforms) & res$Type!='Gene'])
  
  results.mat=rbind(results.mat,c(dataset,types.tab[names(types.tab)=='retained_intron']/sum(types.tab)))
  
  
}


df=data.frame(fraction=as.numeric(as.character(results.mat[1:16,2])),dataset=results.mat[1:16,1],type=ifelse(1:16 %in% c(1:7,15:16) ,'Coronavirus','Other'))

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

wilcox.test(df$frac[c(1:7,15,16)],df$frac[8:14],paired=FALSE)

