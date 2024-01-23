library(glmnet) 
library(rms) 
library(VIM) 
library(survival) 
#读取数据集
df = read.csv("gene_cox.csv")
df=as.matrix(df)
rownames(df)=df[,1]
exp=df[,2:ncol(df)]
dimnames=list(rownames(exp),colnames(exp))
td=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
dt=as.data.frame(td)
for(i in names(dt)[c(3:52)]) {dt[,i] <- as.factor(dt[,i])}
x <- data.matrix(dt[,3:52])
y <- Surv(dt$surtime,dt$surstat)
fit <-glmnet(x,y,family = "cox",alpha = 1)
plot(fit,label=T)
plot(fit,xvar="lambda",label=F,cex.lab=1.2, cex.axis=1, font=1.2, cex.sub=10)#lab 标签 axis 坐标轴 font 加粗
fitcv <- cv.glmnet(x,y,family="cox", alpha=1,nfolds=10)
plot(fitcv,cex.lab=1.2, cex.axis=1, font=1.2)
coef(fitcv, s="lambda.min")


