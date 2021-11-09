library(biomaRt)

source('/Users/karleg/get_fdr_prob.R')

merged.tab=NULL

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host="www.ensembl.org")


for (dataset in c("mason_latest",
                  "SRP040070_3","SRP040070_7","SRP040070_9",
                  "SRP227272_38",'SRP279203_72','SRP294125_74','SRP284977_76','SRP284977_77',"SRP078309_53","SRP178454_50","SRP186406_51","SRP216763_55","SRP222569_54","SRP251704_52",
                  "SRP273785_56"))
{
  res=read.table(paste0('/Users/karleg/HBA-DEALS-covid19_output/',dataset,'.txt'),header=TRUE)
  
  if (dataset!='mason_latest')
  {
    gene.names=getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                     filters    = "ensembl_gene_id",
                     values     = res$Gene, 
                     mart       = mart)
    
    gene.labels=unlist(lapply(res$Gene,function(g)gene.names$hgnc_symbol[gene.names$ensembl_gene_id==g][1]))
    
    res$Gene=gene.labels
  }
  
  res=res[res$Isoform!='Expression',]
  
  res=res[!is.na(res$Gene),]
  
  combined.p=unlist(lapply(res$Gene,function(x)min(res$P[res$Gene==x])))
  
  res$P=combined.p
  
  res=res[!duplicated(res$Gene),]
  
  if (is.null(merged.tab))
  {
    merged.tab=res
  }else{
    merged.tab=merge(merged.tab,res,by='Gene',all=T)
  }
  
}  

merged.tab=merged.tab[,!grepl('Isoform',colnames(merged.tab))]

merged.tab=merged.tab[,!grepl('Explog',colnames(merged.tab))]

ranking=rowSums(1-merged.tab[,2:10],na.rm = T)-rowSums(1-merged.tab[,11:17],na.rm = T)

merged.tab=merged.tab[order(ranking,decreasing = T),]

output=cbind(merged.tab$Gene,ranking[order(ranking,decreasing = T)])

colnames(output)=c('Gene','Score')

write.table(output,'Supp_Table_S2.txt',sep='\t',col.names = T,row.names = F,quote = F)

