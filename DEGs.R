rm(list=ls())
library(limma)
library(edgeR)
df<-read.csv(file="gene_expression.csv")
df=as.matrix(df)
rownames(df)=df[,1]
exp=df[,2:ncol(df)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
group<-c(rep("control",50),rep("treat",374))
design<-model.matrix(~0+factor(group))
colnames(design)<-c("control","treat")
rownames(design)<-colnames(data)
DGElist<-DGEList(counts=data,group=group)
DGElist<-calcNormFactors(DGElist)
v<-voom(DGElist,design,plot=TRUE,normalize="quantile")
contrast.matrix<-makeContrasts(control-treat,levels = design)
fit<-lmFit(v,design)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
DEG<-topTable(fit2,coef="control - treat",n=Inf,adjust="BH")

library(dplyr)
DEG %>%
  mutate(change=case_when(
    logFC >= 1.5 & adj.P.Val<= 0.05 ~"up",
    logFC <= -1.5 & adj.P.Val<= 0.05 ~"down",
    TRUE~"NOT-CHANGE"
  )) -> dif
write.csv(dif,file="DEGs",quote=F)


