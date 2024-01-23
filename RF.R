rm(list=ls())
library(randomForest)
dt=read.csv("gene_cox.csv") 
rownames(dt)=dt[,1]
exp=dt[,2:ncol(dt)]
dt=exp
dt[,3:ncol(dt)]=log2(dt[,3:ncol(dt)]+1)
dt$state <- ifelse(dt$surstat==1 ,1,0)
require(caret)
trainset=dt
folds <- createFolds(y=trainset$state,k=10)
for(i in 1:10){
  fold_test <- trainset[folds[[i]],]   #取folds[[i]]作为测试集
  fold_train <- trainset[-folds[[i]],]   # 剩下的数据作为训练集
  fold_pre <- randomForest(state ~ HGF+TREM2+NQO1+HMOX1+PPARGC1A+EZH2+G6PD+MAPT+MSRA+CDK1+ACADS+ACADL+ESR1+
                             EHHADH+EPO+CDKN2A+BMP6+EGF+CDKN3+AR+RAD51+MMP1+E2F1+S100A8+IGF2BP1+TTR+BIRC5+GSTZ1+
                             CCNA2+CDC25C+TOP2A+MYCN+FOXM1+CCNB1+CPS1+OTC+FANCI+IGF2BP3+PDE2A+PHGDH+ACSL1+SFN+
                             ABAT+MCM2+DMGDH+AURKA+TK1+SARDH+MCM6+FLT3,
                           data = fold_train, importance = TRUE,ntree = 500,mtry=4)
}
fold_predict <- predict(fold_pre,type='response',newdata=fold_test)
fold_predict =ifelse(fold_predict>0.5,1,0)
fold_test$predict = fold_predict
df_fold_test <- data.frame(fold_test)
pred <-prediction(df_fold_test$predict,df_fold_test$state)
perf <- performance(pred,"tpr","fpr")
plot(perf, col='blue',lty=20)
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("ROC curve (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
importance_otu.scale <- data.frame(importance(fold_pre, scale = TRUE), check.names = FALSE)
importance_otu.scale
importance_otu.scale <- importance_otu.scale[order(importance_otu.scale$`%IncMSE`, decreasing = TRUE), ]
importance_otu.scale
library(ggplot2)
imp <- cbind.data.frame(Feature=rownames(fold_pre$importance),fold_pre$importance)
g <- ggplot(imp, aes(x=reorder(Feature, IncNodePurity), y=IncNodePurity))
g + geom_bar(stat = 'identity',fill="lightblue") + xlab('Feature') + theme(axis.text.y = element_text(size = 30)) + theme(axis.title = element_text(size = 30))





