library(limma)
library(VIM)
df=read.csv("gene.csv") 
df=df[,1:7]
df=as.matrix(df)
rownames(df)=df[,1]
exp=df[,2:ncol(df)]
dimnames=list(rownames(exp),colnames(exp))
td=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
td=avereps(td)
td=na.omit(td)
td=as.data.frame(td)
library(survival)
pFilter=0.05 
outResult=data.frame() 
sigGenes=c("surtime","surstat")
for(i in colnames(td[,3:ncol(td)])){ 
  print(i)
  s=Surv(surtime, surstat) ~ td[,i]
  tdcox <- coxph(s, data = td)
  tdcoxSummary = summary(tdcox) 
  pvalue=tdcoxSummary$coefficients[,5] 
  print(pvalue)
  sigGenes=c(sigGenes,i)
  outResult=rbind(outResult,
                  cbind(id=i,
                        HR=tdcoxSummary$conf.int[,"exp(coef)"],
                        L95CI=tdcoxSummary$conf.int[,"lower .95"],
                        H95CI=tdcoxSummary$conf.int[,"upper .95"],
                        pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])
                  )
}
write.csv(outResult,file="临床数据unicox结果1.csv")
UniCoxSurSigGeneExp=td[,sigGenes]
UniCoxSurSigGeneExp=cbind(id=row.names(UniCoxSurSigGeneExp),UniCoxSurSigGeneExp)#以id也就是样品名命名行名
write.csv(UniCoxSurSigGeneExp,file="临床数据unicox结果2.csv")

