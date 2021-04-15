source("git/Entropy/ML_Functions.R")
source("git/Entropy/Preprocess_Functions.R")

BMR_Train = ML.exec(train)
BMR_FCBF = ML.exec(pp)
BMR_KW = ML.exec(Genus_train_KW)
BMR_LDM = ML.exec(Genus_train_LDM)

library(glmnet)
library(caret)
cols = sapply(Genus_train_FCBF, is.numeric)
dataf= Genus_train_FCBF[cols]
targets = as.factor(Genus_train_FCBF$target)
retain = colnames(dataf)


test = readRDS("projects/Entropy/data/test/BV_Genus_AR_test.rds")
ttest = test$target
#makemodel
fit = glmnet(data.matrix(dataf), targets, family = "multinomial", type.multinomial = "grouped")
plot(fit, xvar = "lambda", label = TRUE, type.coef = "2norm")

cvfit=cv.glmnet(data.matrix(dataf), targets, family="multinomial",type.measure = "class", type.multinomial = "grouped", parallel = TRUE,nfolds = 3)
plot(cvfit)

a = predict(cvfit, data.matrix(dataz), s = "lambda.min", type = "class")

cm = confusion.glmnet(cvfit, newx = data.matrix(dataz), newy = ttest, family = "multinomial", s ="lambda.min")


