#!/bin/env Rscript

library(DESeq2)
library(ggplot2)

datalist=paste(c("d11-scr3-1","d11-scr3-2","d11-BACH2-sh2-1","d11-BACH2-sh2-2"),"count.bed",sep='.')

names(datalist)=gsub(".count.bed","",datalist)
data=read.table(datalist[1],sep='\t',stringsAsFactors=F,header=T)
rownames(data)=data[,1]
#data=data.frame()

for(f in 2:length(datalist)){

df1=read.table(datalist[f],sep='\t',stringsAsFactors=F,header=T)
rownames(df1)=df1[,1]
df1=df1[rownames(data),]
data=cbind(data,df1[,2])
}


colnames(data)[2:ncol(data)]=names(datalist)

data$Geneid=NULL
df=data
df=df[which(rowSums(df)!=0),]

# select protein coding gene only 


metadata=data.frame(id=1:ncol(df),dex=rep(c("Control","KO"),each=2),celltype='',geo_id='')
metadata$dex=as.factor(metadata$dex)
#keep=rowSums(df)>=10

#countdata=cbind(geneid=rownames(df[keep,]),df[keep,])
countdata=cbind(geneid=rownames(df),df)

dds <- DESeqDataSetFromMatrix(countData=countdata,
                              colData=metadata,
                              design=~dex, tidy = TRUE)
#use deseq function
dds <- DESeq(dds)

# this will give you the results, P, logFC
res=results(dds)
genelist=c("BACH2","HBB","HBG1","HBG2","HBA1","HBA2")
res[genelist,]
# got the normalized read
norm=data.frame(counts(dds, normalized=TRUE))
res$log10p=-log10(res$pvalue)

threshold=0.05
diff=res[which(res$pvalue<threshold),]
diff$type=''
diff[which(diff$log2FoldChange>0.5),'type']='up'
diff[which(diff$log2FoldChange<(-0.5)),'type']='down'
diff=diff[which(diff$type!=''),]
diff=data.frame(diff)

res$Type='unchange'
res[rownames(diff),'Type']=diff$type

norm=norm[rownames(df),]
res=res[rownames(df),]
data=cbind(df,norm,res)


outfile=paste('DEGs.p=',i,'.2fold.scatter.png',sep='')
png(outfile,width=7*500,height=7*500,res=500)
#anno=data.frame(x=c(-10,10),y=c(40,40),label=c('DOWN: 568','UP: 620'))
#anno=data.frame(x=c(-10,10),y=c(40,40),label=c('DOWN: 81','UP: 319'))

p=ggplot()+
geom_point(res,mapping=aes(log2FoldChange,log10p),color='grey60',size=.5)+
geom_point(diff,mapping=aes(log2FoldChange,log10p,color=type),size=.5)+
xlab("Log2FC")+ylab("-Log10(pvalue)")+
ggtitle("BACH2-sh vs Control")+
#scale_color_gradientn(colours = col,breaks=fill_range,space = "Lab")+
scale_color_manual(values=c('blue','red'))+
coord_cartesian(ylim = c(0, 10))+
#geom_text(data=anno,aes(x,y,label=label),size=7,color="black")+
theme_bw()+
theme(  legend.title = element_blank(),
        legend.text = element_text(size=12,face="bold"),
        legend.position = 'none') +
theme(  axis.title=element_text(color='black', size=18,face="bold"),
        axis.text = element_text(color='black', size=10,face="bold"),
        axis.text.x = element_text(angle = 0,size=9,face='bold'),
        axis.text.y=element_text(angle=0,size=9,face='bold'),
        plot.title = element_text(size=16, face='bold'))+
theme(axis.line = element_line(colour = "black"),
        panel.border=element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold"))+
#geom_abline(slope = k,intercept=0,color="grey60",size=1)+
#geom_vline(xintercept=0,color='grey60',linetype='dashed')+
geom_vline(xintercept=.5,color='red',linetype='dashed',size=.8)+
geom_vline(xintercept=-0.5,color='blue',linetype='dashed',size=.8)+
geom_hline(yintercept=-log10(i),color='grey60',linetype='dashed',size=.8)
print(p);dev.off()

png("genelist.png",,width=7*500,height=7*500,res=500)
anno=data.frame(x=res[genelist,'log2FoldChange'],y=res[genelist,'log10p'],label=genelist)
p=ggplot()+
geom_point(res,mapping=aes(log2FoldChange,log10p),color='grey60',size=.5)+
geom_point(diff[genelist,],mapping=aes(log2FoldChange,log10p,color=type),size=.5)+
xlab("Log2FC")+ylab("-Log10(pvalue)")+
ggtitle("BACH2-sh vs Control")+
#scale_color_gradientn(colours = col,breaks=fill_range,space = "Lab")+
scale_color_manual(values=c('blue','red'))+
coord_cartesian(ylim = c(0, 10))+
geom_text(data=anno,aes(x,y,label=label),size=3,color="black")+
theme_bw()+
theme(  legend.title = element_blank(),
        legend.text = element_text(size=12,face="bold"),
        legend.position = 'none') +
theme(  axis.title=element_text(color='black', size=18,face="bold"),
        axis.text = element_text(color='black', size=10,face="bold"),
        axis.text.x = element_text(angle = 0,size=9,face='bold'),
        axis.text.y=element_text(angle=0,size=9,face='bold'),
        plot.title = element_text(size=16, face='bold'))+
theme(axis.line = element_line(colour = "black"),
        panel.border=element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold"))+
#geom_abline(slope = k,intercept=0,color="grey60",size=1)+
#geom_vline(xintercept=0,color='grey60',linetype='dashed')+
geom_vline(xintercept=0.5,color='red',linetype='dashed',size=.8)+
geom_vline(xintercept=-0.5,color='blue',linetype='dashed',size=.8)+
geom_hline(yintercept=-log10(i),color='grey60',linetype='dashed',size=.8)
print(p);dev.off()



















