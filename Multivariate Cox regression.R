rm(list=ls())
library(survival)
td=read.table("gene_expression.txt",header=T, row.names=1)      
tdmultiCox=coxph(Surv(surtime, surstat) ~ "gene_name", data = td)）
tdmultiCoxSum=summary(tdmultiCox)
outResult=data.frame()
outResult=cbind(
  HR=tdmultiCoxSum$conf.int[,"exp(coef)"],
  L95CI=tdmultiCoxSum$conf.int[,"lower .95"],
  H95CIH=tdmultiCoxSum$conf.int[,"upper .95"],
  pvalue=tdmultiCoxSum$coefficients[,"Pr(>|z|)"])
outResult=cbind(id=row.names(outResult),outResult)
write.table(outResult,file="multiCoxClinical.txt",sep="\t",row.names=F,quote=F)
