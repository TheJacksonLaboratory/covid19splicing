# This script creates a 2-dimensional plot of the GO enrichment scores for two GO terms and the sets of DS genes in the all the datasets (Figure 2).
# Its inputs are the GO enrichment results obtained using the script go_enrichment.R


library(ggpubr)


plot.mat=matrix(ncol=3,nrow=0)

colnames(plot.mat)=c('dataset','x','y')

go.1='GO:0006412'

go.2='GO:0003723'

for (dataset in c("mason_latest",
                  "SRP040070_3","SRP040070_7","SRP040070_9",
                  "SRP227272_38",'SRP279203_72','SRP294125_74',"SRP078309_53","SRP178454_50","SRP186406_51","SRP216763_55","SRP222569_54","SRP251704_52",
                  "SRP273785_56",'SRP284977_76','SRP284977_77'))
{
  
  
  ont.tab=read.table(paste0('GO enrichment covid19/table-',dataset,'_ds-Term-For-Term-Benjamini-Hochberg.txt'),header=TRUE)

  x=ont.tab$Study.term[ont.tab$ID==go.1]
  
  x=x/ont.tab$Study.total[ont.tab$ID==go.1]
  
  x=x*ont.tab$Pop.total[ont.tab$ID==go.1]
  
  x=x/ont.tab$Pop.term[ont.tab$ID==go.1]
  
  x=-log10(ont.tab$p.adjusted[ont.tab$ID==go.1])
  
  if (!go.1 %in% ont.tab$ID)
  {
    x=0
  }else{
    xlab=ont.tab$name[ont.tab$ID==go.1]
  }
  y=ont.tab$Study.term[ont.tab$ID==go.2]
  
  y=y/ont.tab$Study.total[ont.tab$ID==go.2]
  
  y=y*ont.tab$Pop.total[ont.tab$ID==go.2]
  
  y=y/ont.tab$Pop.term[ont.tab$ID==go.2]
  
  y=-log10(ont.tab$p.adjusted[ont.tab$ID==go.2])
  
  if (!go.2 %in% ont.tab$ID)
  {
    y=0
  }else{
    ylab=ont.tab$name[ont.tab$ID==go.2]
  }
  plot.mat=rbind(plot.mat,c(dataset,x,y))
  
}  

betacoronas=c("mason_latest","SRP040070_3","SRP040070_7","SRP040070_9","SRP227272_38",'SRP279203_72','SRP294125_74','SRP284977_76','SRP284977_77')

plot(plot.mat[,2],plot.mat[,3],col=ifelse(plot.mat[,1] %in% betacoronas,'red','blue'),xlab=xlab,ylab=ylab)

labels=c('SARSCOV2','SARS (24h)','MERS (24h)','MERS (48h)','MERS (24h)','SARSCOV2','SARSCOV2','DENV',
         'Strep','HCV','H3N2','H3N2','Zika','RSV','NSP1','NSP2')

df=data.frame(g1=as.numeric(plot.mat[,2]),g2=as.numeric(plot.mat[,3]),is.corona=plot.mat[,1] %in% betacoronas,labels=labels)

ggplot(df, aes(x = g1, y = g2,col=rep('white,3'))) + 
  geom_point() + 
  geom_text(data = df, aes(label = labels,col=is.corona),size=3.5) + theme(legend.position="none")+
  theme(panel.border = element_blank(),
        axis.line    = element_line(color='black'))+theme(panel.background = element_blank())+
  scale_color_manual(values=c("blue","red",'white'))+ labs(x =xlab,y=ylab)+
  geom_hline(yintercept = -log10(0.05),col='black',linetype='dotted')+
  geom_vline(xintercept = -log10(0.05),col='black',linetype='dotted')


