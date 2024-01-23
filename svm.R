rm(list=ls())
library(randomForest)
library(ROCR)
dt=read.csv("gene_cox.csv") 
rownames(dt)=dt[,1]
exp=dt[,2:ncol(dt)]
dt=exp
dt[,3:ncol(dt)]=log2(dt[,3:ncol(dt)]+1)
dt$state<- ifelse(dt$surstat==1,"yes","no")
dt$state_num<- ifelse(dt$surstat==1,1,0)
dt$state<-as.factor(dt$state)
require(caret)
trainset=dt
folds <- createFolds(y=trainset$state,k=10)
for(i in 1:10){
  fold_test <- trainset[folds[[i]],]  
  fold_train <- trainset[-folds[[i]],]   
  fitControl <- trainControl(method="none",classProbs = TRUE)
  fold_pre<- train(state ~ HGF+TREM2+NQO1+HMOX1+PPARGC1A+EZH2+G6PD+MAPT+MSRA+CDK1+ACADS+ACADL+ESR1+
                     EHHADH+EPO+CDKN2A+BMP6+EGF+CDKN3+AR+RAD51+MMP1+E2F1+S100A8+IGF2BP1+TTR+BIRC5+GSTZ1+
                     CCNA2+CDC25C+TOP2A+MYCN+FOXM1+CCNB1+CPS1+OTC+FANCI+IGF2BP3+PDE2A+PHGDH+ACSL1+SFN+
                     ABAT+MCM2+DMGDH+AURKA+TK1+SARDH+MCM6+FLT3,
                   data=fold_train,
                   method="svmRadial",
                   trControl = fitControl,
                   verbose=FALSE,
                   tuneGrid=data.frame(sigma = 0.02,C = 1),
                   metric="ROC")
}
fold_predict <- predict(fold_pre,type='raw',newdata=fold_test)
fold_predict =ifelse(fold_predict=="yes",1,0)
fold_test$predict = fold_predict
df_fold_test <- data.frame(fold_test)
pred <-prediction(df_fold_test$predict,df_fold_test$state_num)
perf <- performance(pred,"tpr","fpr")
plot(perf, col='blue',lty=20)
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("ROC curve (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
plot(fold_pre)
fold_predict <- predict(fold_pre,type='raw',newdata=fold_train)
a=varImp(fold_pre)  
a$importance
library(ggplot2)
imp <- cbind.data.frame(Feature=rownames(a$importance),a$importance)
g <- ggplot(imp, aes(x=reorder(Feature, yes), y=yes))
g + geom_bar(stat = 'identity',fill="lightblue") + xlab('Feature') + ylab('Importance')+ theme(axis.text.y = element_text(size = 30)) + theme(axis.title = element_text(size = 30))


