#This script generates the file 'rbp_enrichment_sarscov2.txt', which is then used as input for the script 'analyze_rbp_enrichment.R'
#that generates figure 7.
#The input files for this scripts are:
# 1. HBA-DEALS output from the Snakemake pipeline
# 2. Files containing the ORNAment database

library(data.table)

library(edgeR)

library(biomaRt)

source('get_fdr_prob.R')

enrich.mat=matrix(ncol=4,nrow=0)

colnames(enrich.mat)=c('RBP','Dataset','FC','P.Value')

min.score=0.95

encoding.tab=read.csv('Homo_sapiens_string_to_int_ID_conversion.csv',header = F)

colnames(encoding.tab)=c('ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name', 'ensembl_transcript_id_INT', 'ensembl_gene_id_INT')

skip=0

ornament=matrix(ncol=12,nrow=0)

next.read=matrix(ncol=12,nrow=0)

colnames(ornament)=c('ensembl_gene_id', 'ensembl_transcript_id', 'gene_biotype', 'transcript_biotype', 'transcript_position', 'RBP', 'score', 'unpaired_probability', 
                     'chromosome', 'region', 'exon_start', 'exon_end')

prev.rows=-1

while(prev.rows<=nrow(next.read))
{
  
  prev.rows=nrow(next.read)
  
  next.read=fread('Homo_sapiens_cDNA_oRNAment.csv',nrows = 10000000,header=FALSE,skip=skip,data.table = FALSE,nThread = 40)
  
  ornament=rbind(ornament,as.matrix(next.read))
  
  ornament=ornament[ornament[,'score']>=min.score,]
  
  ornament=ornament[!duplicated(paste(ornament[,2],ornament[,6])),]
  
  skip=skip+10000000
  
}

rbp.tab=read.csv('RBP_id_encoding.csv',header=FALSE)

colnames(rbp.tab)=c('rbp.num','RBP')

rbp.tab$RBP=unlist(lapply(unlist(lapply(rbp.tab$RBP,function(x)gsub(' ','',x))),function(x)gsub('\\(.*\\)','',x)))

ornament=merge(ornament,rbp.tab,by.x='RBP',by.y='rbp.num')

ornament$ensembl_transcript_id=as.integer(ornament$ensembl_transcript_id)

ornament=merge(ornament,encoding.tab,by.x='ensembl_transcript_id',by.y='ensembl_transcript_id_INT')

for (dataset in c(
                  "SRP040070_3","SRP040070_7","SRP040070_9",
                  "SRP227272_38",'SRP279203_72','SRP294125_74',"SRP078309_53","SRP178454_50","SRP186406_51","SRP216763_55","SRP222569_54","SRP251704_52",
                  "SRP273785_56",'SRP284977_76','SRP284977_77'))
{   
  res=read.table(paste0('HBA-DEALS-covid19_output/',dataset,'.txt'),header=TRUE)
  
  thresh=get.fdr.prob(res)
  
  ds.isoforms=res$Isoform[res$Isoform!='Expression' & res$P<=thresh[[2]]]
  
  next.ornament=ornament[ornament$ensembl_transcript_id.y %in% res$Isoform,]
  
  for (rbp.id in unique(next.ornament$RBP))
  {
    
    white.balls.drawn=sum(next.ornament$RBP==rbp.id & (next.ornament$ensembl_transcript_id.y %in% ds.isoforms))
    
    white.balls.present=sum(unique(next.ornament$ensembl_transcript_id.y) %in% ds.isoforms)
    
    black.balls.present=length(unique(next.ornament$ensembl_transcript_id.y))-white.balls.present
    
    number.draws=sum(next.ornament$RBP==rbp.id)
    
    fc=1
    
    if (number.draws>0 && white.balls.present>0)
    
        fc=white.balls.drawn/number.draws*(white.balls.present+black.balls.present)/white.balls.present
    
    if (fc>1)
    {
        p=phyper(q = white.balls.drawn-1,m = white.balls.present,n=black.balls.present,k = number.draws,lower.tail = F )
    }else{
      
        p=phyper(q = white.balls.drawn,m = white.balls.present,n=black.balls.present,k = number.draws)
    }
    
    enrich.mat=rbind(enrich.mat,c(rbp.id,dataset,fc,p))
    
  }
  write.table(enrich.mat,'rbp_enrichment_sarscov2.txt',sep='\t',col.names = T,row.names = F,quote = F)
  
}
