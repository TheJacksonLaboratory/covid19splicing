# This script runs Ontologizer on sets of differentially spliced and differentially expressed genes in order to find enriched
# GO categories.  The .obo and .gaf files that Ontologizer uses as inputs can be obtained from the GO website.  The output if this
# script is stored in the directory 'GO enrichment covid19' for use by downstream scripts.

library(biomaRt)

source('get_fdr_prob.R')

for (dataset in c("mason_latest",
                  "SRP040070_3","SRP040070_7","SRP040070_9",
                  "SRP227272_38",'SRP279203_72','SRP294125_74',"SRP078309_53","SRP178454_50","SRP186406_51","SRP216763_55","SRP222569_54","SRP251704_52",
                  "SRP273785_56",'SRP284977_76','SRP284977_77'))
{
  res=read.table(paste0('HBA-DEALS-covid19_output/',dataset,'.txt'),header=TRUE)
  
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host="www.ensembl.org")
  
  gene.names=getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                     filters    = "ensembl_gene_id",
                     values     = res$Gene, 
                     mart       = mart)
  
  gene.labels=unlist(lapply(res$Gene,function(g)gene.names$hgnc_symbol[gene.names$ensembl_gene_id==g][1]))
  
  if (dataset=='mason_latest')
    
    gene.labels=res$Gene
  
  thresh=get.fdr.prob(res)
  
  de.genes=gene.labels[res$Isoform=='Expression' & res$P<=thresh[[1]]]
  
  ds.genes=unique(gene.labels[res$Isoform!='Expression' & res$P<=thresh[[2]]])
  
  dast.genes=ds.genes
  
  dge.genes=de.genes
  
  all.genes<-unique(gene.labels)
  
  write.table(all.genes,'universe.txt',quote = F,row.names = F,col.names = F)
  
  write.table(dge.genes,paste0(dataset,'_de.txt'),quote = F,row.names = F,col.names = F)
  
  system(paste0('java -jar Ontologizer.jar -g go.obo -a goa_human.gaf -s ',dataset,'_de.txt -p universe.txt -c Term-For-Term -m Benjamini-Hochberg -n'))
  
  
  write.table(dast.genes,paste0(dataset,'_ds.txt'),quote = F,row.names = F,col.names = F)
  
  system(paste0('java -jar Ontologizer.jar -g go.obo -a goa_human.gaf -s ',dataset,'_ds.txt -p universe.txt -c  Term-For-Term -m Benjamini-Hochberg -n'))
}  
