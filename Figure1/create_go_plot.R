
# This script uses the GO enrichment results fron the script 'go_enrichment.R' to plot figure 1.  Its input also includes the file
# featured_GO_terms.tsv which contains the GO terms to plot and is located in the same directory as this script.

library(ggplot2)
library(gridExtra)
library(ggsci)

datasets=c("mason_latest","SRP040070_3","SRP040070_7","SRP040070_9",
                           "SRP227272_38",'SRP279203_72','SRP294125_74','SRP284977_77','SRP284977_76',"SRP078309_53","SRP178454_50","SRP186406_51","SRP216763_55","SRP222569_54","SRP251704_52",
           "SRP273785_56") 

names.all=c('SARS-COV2(A)','SARS-High-24h','MERS-High-24h','MERS-High-48h','MERS','SARS-COV2(B)',"SARS-COV2(C)",'NSP2','NSP1','Dengue','S.Pneumoniae','HCV','H3N2',"H3N2'",'Zika','RSV')

GO_gp=as.data.frame(matrix(nrow=0,ncol=5,dimnames = list(c(),c("Group","GO_biological_process","Gene_number","Fold_enrichment","FDR" ))))

for (dataset in datasets)
{
  Ontologizer.table=read.table(paste0('GO enrichment covid19/table-',dataset,'_ds-Term-For-Term-Benjamini-Hochberg.txt'),header=TRUE)
 
  if (sum(Ontologizer.table$p.adjusted<=0.05)==0)
    
    next
  
  GO_all=data.frame(GO_biological_process=paste0(Ontologizer.table$name,'(',Ontologizer.table$ID,')'))
  
  GO_all=cbind(GO_all,Ontologizer.table$Study.term)
  
  GO_all=cbind(GO_all,Ontologizer.table$Study.term/Ontologizer.table$Study.total/Ontologizer.table$Pop.term*Ontologizer.table$Pop.total)
  
  GO_all=cbind(GO_all,Ontologizer.table$p.adjusted)
  
  colnames(GO_all)=c("GO_biological_process","Gene_number","Fold_enrichment","FDR" )
  
  GO_all=cbind(rep(names.all[which(datasets==dataset)],nrow(GO_all)),GO_all)
  
  colnames(GO_all)[1]='Group'
  
  GO_gp=rbind(GO_gp,GO_all)
  
  
}

feature.go.terms=read.csv('featured_GO_terms.tsv',header = T,sep='\t')

GO_gp=GO_gp[GO_gp$GO_biological_process %in% paste0(feature.go.terms$GO.label,'(',feature.go.terms$GO.id,')'),]

GO_gp$GO_biological_process[GO_gp$GO_biological_process=='protein−containing complex assembly(GO:0065003)']='PCCA(GO:0065003)'

GO_gp$GO_biological_process[GO_gp$GO_biological_process=='protein localization to endoplasmic reticulum(GO:0070972)']='PLTER(GO:0070972)'

GO_gp$GO_biological_process[GO_gp$GO_biological_process=='biological process involved in interspecies interaction between organisms(GO:0044419)']=
  'BPIIIIBO(GO:0044419)'

GO_gp$GO_biological_process[GO_gp$GO_biological_process=='cellular response to stress(GO:0033554)']=
  'CRTS(GO:0033554)'

GO_gp$GO_biological_process[GO_gp$GO_biological_process=='posttranscriptional regulation of gene expression(GO:0010608)']='PROGE(GO:0010608)'

GO_gp$GO_biological_process[GO_gp$GO_biological_process=='negative regulation of biosynthetic process(GO:0009890)']='NROBP(GO:0009890)'

# List objects and their structure contained in the dataframe 'GO_gp'
ls.str(GO_gp)

# Transform the column 'Gene_number' into a numeric variable
GO_gp$Gene_number <- as.numeric(GO_gp$Gene_number)

# Replace all the "_" by a space in the column containing the GO terms
GO_gp$GO_biological_process <- chartr("_", " ", GO_gp$GO_biological_process)

# Transform the column 'GO_biological_process' into factors
GO_gp$GO_biological_process<-as.factor(GO_gp$GO_biological_process)

# Transform FDR values by -log10('FDR values')
GO_gp$'|log10(FDR)|' <- -(log10(GO_gp$FDR))

# Change factor order
GO_gp$Group<- factor(GO_gp$Group,levels = names.all)
GO_gp$GO_biological_process<-factor(GO_gp$GO_biological_process,levels=rev(levels(GO_gp$GO_biological_process)))

# Create a vector with new names for groups to use in the plot
# Replace the terms by your own (\n allow to start a new line)
group.labs <- c(`SARS-COV2(A)`="SARS-COV2(A)",`SARS-High-24h` = "SARS-High-24h",
                `MERS-High-24h` = "MERS-High-24h",
                `MERS-High-48h` = "MERS-High-48h",
                `MERS` = "MERS",
                `SARS-COV2(B)`="SARS-COV2(B)",
                `SARS-COV2(C)`="SARS-COV2(C)",
                  `NSP1`="NSP1",
                  `NSP2`="NSP2",
                `Dengue`="Dengue",
                 `S.Pneumoniae`="S.Pneumoniae",
                 `HCV`="HCV",
                 `H3N2`="H3N2",
                 `H3N2'`="H3N2'",
                 `Zika`="Zika",
                 `RSV`="RSV"
                )

GO_gp$GO_biological_process=gsub('\\(.*\\)','',GO_gp$GO_biological_process)

ggplot(GO_gp, aes(x = GO_biological_process, y = Fold_enrichment)) +
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "azure4", size=.5)+
  geom_point(data=GO_gp,aes(x=GO_biological_process, y=Fold_enrichment,size = Gene_number, colour = `|log10(FDR)|`), alpha=.7)+
  scale_color_gradient(low="blue",high="red",limits=c(0, NA))+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt")),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
        axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        legend.title.align=0.5)+
  xlab("GO term")+
  ylab("Fold enrichment")+
  labs(color="-log10(FDR)", size="Number\nof genes")+
  facet_wrap(~Group,ncol=4,labeller=as_labeller(group.labs))+
  theme(text = element_text(size=15))+scale_fill_npg()+ theme(axis.text = element_text(size = 14))+ theme(axis.title.x  = element_text(size = 16)) 

