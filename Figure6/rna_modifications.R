#This script generates figure 6.  Its input is the RMBase database files and the HBA-DEALS results from the Snakemake pipeline.

source('get_fdr_prob.R')

library(biomaRt)

library(ggpubr)

res.mat=matrix(ncol=4,nrow=0)

colnames(res.mat)=c('Database','modifications','fc','pval')

mod.types=c('m6A_site','PseudoU_site','otherMod_site','m5C_site','Nm_site','mod_RBP_eraser','mod_RBP_writer','m1A_site','mod_RBP_reader')

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host="www.ensembl.org")

for (mod.type in mod.types)
{
  
  next.mod=read.table(paste0('RMBase_hg19_all_',mod.type,'.txt'),sep='\t')
  
  for (dataset in c("mason_latest",
    "SRP040070_3","SRP040070_7","SRP040070_9",
    "SRP227272_38",'SRP279203_72','SRP294125_74',"SRP078309_53","SRP178454_50","SRP186406_51","SRP216763_55","SRP222569_54","SRP251704_52",
    "SRP273785_56",'SRP284977_76','SRP284977_77'))
  {   
    res=read.table(paste0('HBA-DEALS-covid19_output/',dataset,'.txt'),header=TRUE)
    
    thresh=get.fdr.prob(res)
    
    gene.names=getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                     filters    = "ensembl_gene_id",
                     values     = res$Gene, 
                     mart       = mart)
    
    gene.labels=unlist(lapply(res$Gene,function(g)gene.names$hgnc_symbol[gene.names$ensembl_gene_id==g][1]))
    
    if (dataset=='mason_latest')
      
      gene.labels=res$Gene
    
    ds.genes=unique(gene.labels[res$Isoform!='Expression' & res$P<=thresh[[2]]])
   
    white.balls.drawn=sum(ds.genes %in% unlist(lapply(next.mod$V12,strsplit,',')))
    
    white.balls.present=length(ds.genes)
    
    black.balls.present=length(unique(res$Gene))-white.balls.present
    
    number.draws=sum(unique(gene.labels) %in% unlist(lapply(next.mod$V12,strsplit,',')))
    
    fc=1
    
    if (number.draws>0 && white.balls.present>0)
      
      fc=white.balls.drawn/number.draws*(white.balls.present+black.balls.present)/white.balls.present
    
    if (fc>1)
    {
      p=phyper(q = white.balls.drawn-1,m = white.balls.present,n=black.balls.present,k = number.draws,lower.tail = F )
    }else{
      
      p=phyper(q = white.balls.drawn,m = white.balls.present,n=black.balls.present,k = number.draws)
    }
    
    
    res.mat=rbind(res.mat,c(dataset,mod.type,fc,p))
    
     
  }
}

betacoronas=c("mason_latest","SRP040070_3","SRP040070_7","SRP040070_9","SRP227272_38",'SRP279203_72','SRP294125_74','SRP284977_76','SRP284977_77')

df=data.frame(modification=res.mat[,2],is.betacorona=res.mat[,1] %in% betacoronas,minus.log10.pval=-log10(as.numeric(res.mat[,4])))

ggboxplot(data = df,x = 'modification',y = 'minus.log10.pval',palette = 'npg',fill = 'is.betacorona',merge = TRUE,add='jitter',shape='is.betacorona')+
  geom_hline(yintercept = -log10(0.05),col='red',linetype='dotted')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.text = element_text(size = 20))+ theme(axis.title.y  = element_text(size = 20))+ theme(axis.title.x  = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+theme(legend.title = element_text(size = 20))+ylab('-log10(p-value)')

wilcox.test(df$minus.log10.pval[df$modification=='PseudoU_site' & df$is.betacorona==TRUE],
df$minus.log10.pval[df$modification=='PseudoU_site' & df$is.betacorona==FALSE],paired = FALSE)


df=data.frame(modification=res.mat[,2],is.betacorona=res.mat[,1] %in% betacoronas,FC=as.numeric(res.mat[,3]))

ggboxplot(data = df,x = 'modification',y = 'FC',palette = 'npg',fill = 'is.betacorona',merge = TRUE,add='jitter',shape='is.betacorona')+
  geom_hline(yintercept = 1,col='red',linetype='dotted')+ggtitle('Fold-change')+theme(axis.text.x = element_text(angle = 90, hjust = 1))

