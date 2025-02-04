set.seed(11)
setwd("C:/Users/HayBa/Documents/School/2020-21/Rotation 2/R/MethylationDE/FIXED/")

library(tidyverse)
library(sesame)
library(limma)
library(readr)
library(ggplot2)
library(corrplot)
library(qqman)
library(caret)
library(dplyr)
library(DESeq2)
library(plotROC)
library(ggbeeswarm)
library(rpart)
library(shapper)
library(RVenn)


###### processing beta values ######
betas = read_csv("CADET_rawbetas.csv")
betas = as.data.frame(betas)
rownames(betas) <- betas$X1
betas$X1 <- NULL

# remove problematic probes
missingProbe <- apply(betas, 1, function(x) sum(is.na(x)))
betas.clean <- betas[missingProbe < ncol(betas), ]
nrow(betas.clean) #759,450

# remove probes w/ 25% missingness or greater using mean imputation for remaining NAs
missingProbe <- apply(betas.clean, 1, function(x) sum(is.na(x)))
betas.clean <- betas.clean[missingProbe <= .25 * ncol(betas.clean), ]
nrow(betas.clean) #750,582

# remove non-CpG probes (non-CpG methylation site ("ch") and SNPs("rs"))
betas.clean = betas.clean[-grep("ch", rownames(betas.clean)), ]
betas.clean = betas.clean[-grep("rs", rownames(betas.clean)), ]
nrow(betas.clean) #748,275

# now impute
k = which(is.na(betas.clean), arr.ind = T)
betas.clean[k] = rowMeans(betas.clean, na.rm = T)[k[,1]]

# logit transform beta values to M values (still bimodally distributed, better approximate normal dist.)
plot(density(as.matrix(betas.clean)))
Ms = BetaValueToMValue(betas.clean)
plot(density(as.matrix(Ms)))

######### PCA to identify outliers and patterns in data #######
# load phenotype data
pd <- read.csv("CADET_SampleInfoMethylation_SingleFamily_12_6.csv")

# remove individuals with mising age or sex
table(pd$Age, useNA = "always")
table(pd$Sex, useNA = "always")
pd <- pd[pd$Sex != "", ]
nrow(pd) #172

Ms = Ms[ , colnames(Ms) %in% pd$Basename]
dim(Ms) #748275 172

pd <- pd[pd$Basename %in% colnames(Ms), ]
dim(pd) #172 15

# samples separated by Sex because includes probes on sex chromosomes without adjusting
# remove sex chromosome probes and try again
# read in array manifest
ann <- read.table("EPIC.hg19.manifest.tsv.gz", sep = "\t", header = T)

Ms_nosex <- Ms[rownames(Ms) %in% ann[ann$CpG_chrm != "chrX" & ann$CpG_chrm != "chrY", ]$probeID, ]
dim(Ms_nosex) #732300 172

# only want PrePF
pd <- pd[pd$DX != "IPF", ]
rownames(pd) <- pd$Basename
nrow(pd) #130
#Ms_nosex <- Ms_nosex[ , colnames(Ms_nosex) %in% rownames(pd)]
Ms <- Ms[ , colnames(Ms) %in% rownames(pd)]

dat <- t(Ms)
dat <- merge(pd %>% select(DX), dat, by = "row.names")
rownames(dat) <- dat$Row.names
dat$Row.names <- NULL
# dat is now samples in rows and methylation sites in columns

control <- trainControl(method = "repeatedcv", #k-folds cross-validation
                        number = 5, #5 folds
                        repeats = 10, #repeat 10 times
                        classProbs = T, #class probability predictions for classification
                        summaryFunction = twoClassSummary, #performance metrics for 2 class problems
                        savePredictions = T) #saves predictions for each resample

col_var <- colVars(as.matrix(dat[ ,2:ncol(dat)]))

prepfvals <- read.csv("prepf.conts.cadet.glint.lmm.txt", row.names = 1)
prepfvals <- prepfvals[order(prepfvals$p.values), ]

################ glmnet 1000 features by differential testing ##############
sorted <- sort(col_var, decreasing = T, index.return = T)$ix[1:1000]
sorted <- sorted + 1
dat_sub <- dat[, c(1,sorted)]
colnames(dat_sub) <- make.names(colnames(dat_sub))

#relevel for caret (takes first level as positive case)
dat_sub$DX <- as.factor(dat_sub$DX)
dat_sub$DX <- relevel(dat_sub$DX, ref = "PrePF")

#split data into training and test sets
#split 70:30 while maininting equivalent numbers of each DX in both groups
inTrain <- createDataPartition(dat_sub$DX, p = 0.7, list = F)

dat_train <- dat_sub[inTrain, ] #92
dat_test <- dat_sub[-inTrain, ] #38

prepf <- prepfvals[1:1000, ] # subset top 1000

prepf <- dat[, colnames(dat) %in% rownames(prepf)] #subset cpgs in glint
prepf <- cbind(DX = dat$DX, prepf)
#refactor for caret
prepf$DX <- as.factor(prepf$DX)
prepf$DX <- relevel(prepf$DX, ref = "PrePF")

prepf_train <- prepf[rownames(prepf) %in% rownames(dat_train), ]
prepf_test <- prepf[!rownames(prepf) %in% rownames(dat_train), ]

colnames(prepf_train) <- make.names(colnames(prepf_train))
colnames(prepf_test) <- make.names(colnames(prepf_test))


model <- train(DX ~ ., data = prepf_train, method = "glmnet", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$alpha == final_param$alpha & model$results$lambda == final_param$lambda, ]
model_final #good specificity and sensitivity

# try some tuning
tuning <- expand.grid(alpha = c(0.025, 0.05, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 
                                0.115, 0.12, 0.13, 0.14, 0.15, 0.2), 
                      lambda = c(0.01, 0.015, 0.0175, 0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 
                                 0.026, 0.027, 0.028, 0.029, 0.03, 0.04, 0.05))
model2 <- train(DX ~ ., data = prepf_train, method = "glmnet", trControl = control, 
                tuneGrid = tuning)

final_param2 <- model2$bestTune #final model
model_final2 <- model2$results[model2$results$alpha == final_param2$alpha & model2$results$lambda == final_param2$lambda, ]
model_final2 #improved all

#look at top sites in model
top_vars <- varImp(model2)$importance
top_vars$site <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$Overall, decreasing = T), ]
top_vars <- top_vars[top_vars$Overall != 0, ]
top_vars <- top_vars[1:75, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), Overall), y = Overall)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = Overall), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#look at results
predicted <- predict(model2, newdata = prepf_test)
predicted

test_probs <- predict(model2, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX)

#select training predictions using final hyperparamter values
select_idx <- model2$pred$alpha == final_param2$alpha & 
  model2$pred$lambda == final_param2$lambda
g <- ggplot(model2$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model2$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_glmnet1000diff <- model2
prepf_glmnet1000diff_testResults <- test_results

################ SVMLinear 1000 features by differential testing ##################
model <- train(DX ~ ., data = prepf_train, method = "svmLinear", trControl = control, 
               tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$C == final_param$C, ]
model_final #very good model

#try tuning
tuning <- expand.grid(C = c(0.01, 0.1, 1, 5, 10, 15, 20, 25, 50, 100, 150, 200, 500, 1000))
model2 <- train(DX ~ ., data = prepf_train, method = "svmLinear", trControl = control, 
                tuneGrid = tuning)

final_param2 <- model2$bestTune #final model
model_final2 <- model2$results[model2$results$C == final_param2$C, ]
model_final2 #made it worse

plot(model)

#look at results
predicted <- predict(model, newdata = prepf_test)

test_probs <- predict(model, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX)

#look at top sites in model
top_vars <- varImp(model)$importance
top_vars$site <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$PrePF, decreasing = T), ]
top_vars <- top_vars[top_vars$PrePF != 0, ] #994
top_vars <- top_vars[1:75, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), PrePF), y = PrePF)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = PrePF), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#select training predictions using final hyperparamter values
select_idx <- model$pred$C == final_param$C
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_svml1000diff <- model
prepf_svml1000diff_testResults <- test_results

################ svmRadial 1000 features by differential testing ###################
model <- train(DX ~ ., data = prepf_train, method = "svmRadial", trControl = control, 
               tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$sigma == final_param$sigma & model$results$C == final_param$C, ]
model_final #good overall; not the best sensitivity

# try tuning
tuning <- expand.grid(sigma = c(0.00047, 0.00048, 0.00049, 0.0005, 0.00051, 0.00052, 0.00053, 
                                0.00054, 0.00055, 0.00056, 0.0006), 
                      C = c(0.01, 0.1, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 20, 
                            25, 30, 50))
model2 <- train(DX ~ ., data = prepf_train, method = "svmRadial", trControl = control, 
                tuneGrid = tuning)

final_param2 <- model2$bestTune #final model
model_final2 <- model2$results[model2$results$sigma == final_param2$sigma & model2$results$C == final_param2$C, ]
model_final2 #didn't do much lol rip

plot(model)
plot(model2)

#look at results
predicted <- predict(model, newdata = prepf_test)

test_probs <- predict(model, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX)

#look at top sites in model
top_vars <- varImp(model)$importance
top_vars$site <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$PrePF, decreasing = T), ]
top_vars <- top_vars[top_vars$PrePF != 0, ] #994
top_vars <- top_vars[1:75, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), PrePF), y = PrePF)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = PrePF), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#select training predictions using final hyperparamter values
select_idx <- model$pred$C == final_param$C & model$pred$sigma == final_param$sigma
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_svmr1000diff <- model
prepf_svmr1000diff_testResults <- test_results

################ RF 1000 features by differential testing ###############
model <- train(DX ~ ., data = prepf_train, method = "rf", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$mtry == final_param$mtry, ]
model_final #bad sensitivity; weird

plot(model)

top_vars <- varImp(model)$importance
top_vars$gene <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$Overall, decreasing = T), ]
top_vars <- top_vars[top_vars$Overall != 0, ] #877
top_vars <- top_vars[1:75, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), Overall), y = Overall)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = Overall), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#look at results
predicted <- predict(model, newdata = prepf_test)
predicted

#also want predicted class probabilities 
test_probs <- predict(model, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX) 

#select training predictions using final hyperparamter values
select_idx <- model$pred$mtry == final_param$mtry
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_rf1000diff <- model
prepf_rf1000diff_testResults <- test_results
################ glmnet 100 features by differential testing #############
sorted <- sort(col_var, decreasing = T, index.return = T)$ix[1:100]
sorted <- sorted + 1
dat_sub <- dat[, c(1,sorted)]
colnames(dat_sub) <- make.names(colnames(dat_sub))

#relevel for caret (takes first level as positive case)
dat_sub$DX <- as.factor(dat_sub$DX)
dat_sub$DX <- relevel(dat_sub$DX, ref = "PrePF")

#split data into training and test sets
#split 70:30 while maininting equivalent numbers of each DX in both groups
inTrain <- createDataPartition(dat_sub$DX, p = 0.7, list = F)

dat_train <- dat_sub[inTrain, ] #92
dat_test <- dat_sub[-inTrain, ] #38

prepf <- prepfvals[1:100, ]

prepf <- dat[, colnames(dat) %in% rownames(prepf)] #subset cpgs in glint
prepf <- cbind(DX = dat$DX, prepf)
#refactor for caret
prepf$DX <- as.factor(prepf$DX)
prepf$DX <- relevel(prepf$DX, ref = "PrePF")

prepf_train <- prepf[rownames(prepf) %in% rownames(dat_train), ]
prepf_test <- prepf[!rownames(prepf) %in% rownames(dat_train), ]

colnames(prepf_train) <- make.names(colnames(prepf_train))
colnames(prepf_test) <- make.names(colnames(prepf_test))

model <- train(DX ~ ., data = prepf_train, method = "glmnet", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$alpha == final_param$alpha & model$results$lambda == final_param$lambda, ]
model_final #good specificity and sensitivity

# try some tuning
tuning <- expand.grid(alpha = c(0.05, 0.06, 0.07, 0.08, 0.09, 0.0925, 0.095, 0.0975, 0.1, 0.105, 
                                0.11, 0.115, 0.12, 0.13, 0.14, 0.15, 0.2), 
                      lambda = c(0.01, 0.0125, 0.015, 0.0175, 0.018, 0.019, 0.02, 0.021, 0.022,
                                 0.0225, 0.025, 0.0275, 0.03))
model2 <- train(DX ~ ., data = prepf_train, method = "glmnet", trControl = control, 
                tuneGrid = tuning)

final_param2 <- model2$bestTune #final model
model_final2 <- model2$results[model2$results$alpha == final_param2$alpha & model2$results$lambda == final_param2$lambda, ]
model_final2 #made it worse keep model 1

plot(model)

#look at top sites in model
top_vars <- varImp(model)$importance
top_vars$site <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$Overall, decreasing = T), ]
top_vars <- top_vars[top_vars$Overall != 0, ] #72
#top_vars <- top_vars[1:75, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), Overall), y = Overall)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = Overall), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#look at results
predicted <- predict(model, newdata = prepf_test)
predicted

test_probs <- predict(model, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX)

#select training predictions using final hyperparamter values
select_idx <- model$pred$alpha == final_param$alpha & 
  model$pred$lambda == final_param$lambda
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_glmnet100diff <- model
prepf_glmnet100diff_testResults <- test_results

################ SVMLinear 100 features by differential testing ##################
model <- train(DX ~ ., data = prepf_train, method = "svmLinear", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$C == final_param$C, ]
model_final #good specificity and sensitivity overall

# try tuning
tuning <- expand.grid(C = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 0.95, 0.975, 
                            1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 5, 10, 20))
model2 <- train(DX ~ ., data = prepf_train, method = "svmLinear", trControl = control, 
                tuneGrid = tuning)

final_param2 <- model2$bestTune #final model
model_final2 <- model2$results[model2$results$C == final_param2$C, ]
model_final2 #made it worse

#look at results
predicted <- predict(model, newdata = prepf_test)

test_probs <- predict(model, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX)

#look at top sites in model
top_vars <- varImp(model)$importance
top_vars$site <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$PrePF, decreasing = T), ]
top_vars <- top_vars[top_vars$PrePF != 0, ] #994
top_vars <- top_vars[1:75, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), PrePF), y = PrePF)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = PrePF), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#select training predictions using final hyperparamter values
select_idx <- model$pred$C == final_param$C
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_svml100diff <- model
prepf_svml100diff_testResults <- test_results

################ svmRadial 100 features by differential testing ###################
model <- train(DX ~ ., data = prepf_train, method = "svmRadial", trControl = control, 
               tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$sigma == final_param$sigma & model$results$C == final_param$C, ]
model_final #good specificity andgood sensitivity too

# try tuning
tuning <- expand.grid(sigma = c(0.004, 0.0045, 0.005, 0.0051, 0.0052, 0.0053, 0.0054, 0.0055, 
                                0.0056, 0.0057, 0.0056, 0.00575, 0.006), 
                      C = c(2, 3, 4, 5, 6, 7, 8, 10, 25, 50))
model2 <- train(DX ~ ., data = prepf_train, method = "svmRadial", trControl = control, 
                tuneGrid = tuning)

final_param2 <- model2$bestTune #final model
model_final2 <- model2$results[model2$results$sigma == final_param2$sigma & model2$results$C == final_param2$C, ]
model_final2 #improved all 3

# try tuning again
tuning2 <- expand.grid(sigma = c(0.0035, 0.0036, 0.0037, 0.0038, 0.0039, 0.004, 0.0041, 0.0042, 
                                 0.0043, 0.0044, 0.0045), 
                       C = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 50))
model3 <- train(DX ~ ., data = prepf_train, method = "svmRadial", trControl = control, 
                tuneGrid = tuning2)

final_param3 <- model3$bestTune #final model
model_final3 <- model3$results[model3$results$sigma == final_param3$sigma & model3$results$C == final_param3$C, ]
model_final3 #made it worse; keep model 2

plot(model2)

#look at results
predicted <- predict(model2, newdata = prepf_test)

test_probs <- predict(model2, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX)

#look at top sites in model
top_vars <- varImp(model2)$importance
top_vars$site <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$PrePF, decreasing = T), ]
top_vars <- top_vars[top_vars$PrePF != 0, ] #994
top_vars <- top_vars[1:75, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), PrePF), y = PrePF)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = PrePF), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#select training predictions using final hyperparamter values
select_idx <- model2$pred$C == final_param2$C & model2$pred$sigma == final_param2$sigma
g <- ggplot(model2$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model2$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_svmr100diff <- model2
prepf_svmr100diff_testResults <- test_results

################ RF 100 features by differential testing ###############
model <- train(DX ~ ., data = prepf_train, method = "rf", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$mtry == final_param$mtry, ]
model_final #pretty good; bad sensitivity; still sucky lol

top_vars <- varImp(model)$importance
top_vars$gene <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$Overall, decreasing = T), ]
top_vars <- top_vars[top_vars$Overall != 0, ] #877
top_vars <- top_vars[1:75, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), Overall), y = Overall)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = Overall), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

plot(model)

#look at results
predicted <- predict(model, newdata = prepf_test)
predicted

#also want predicted class probabilities 
test_probs <- predict(model, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX) 

#select training predictions using final hyperparamter values
select_idx <- model$pred$mtry == final_param$mtry
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_rf100diff <- model
prepf_rf100diff_testResults <- test_results

################ glmnet 50 features by differential testing ################
sorted <- sort(col_var, decreasing = T, index.return = T)$ix[1:50]
sorted <- sorted + 1
dat_sub <- dat[, c(1,sorted)]
colnames(dat_sub) <- make.names(colnames(dat_sub))

#relevel for caret (takes first level as positive case)
dat_sub$DX <- as.factor(dat_sub$DX)
dat_sub$DX <- relevel(dat_sub$DX, ref = "PrePF")

#split data into training and test sets
#split 70:30 while maininting equivalent numbers of each DX in both groups
inTrain <- createDataPartition(dat_sub$DX, p = 0.7, list = F)

dat_train <- dat_sub[inTrain, ] #92
dat_test <- dat_sub[-inTrain, ] #38

prepf <- prepfvals[1:50, ] # subset top 100

prepf <- dat[, colnames(dat) %in% rownames(prepf)] #all three
prepf <- cbind(DX = dat$DX, prepf)
#refactor for caret
prepf$DX <- as.factor(prepf$DX)
prepf$DX <- relevel(prepf$DX, ref = "PrePF")

prepf_train <- prepf[rownames(prepf) %in% rownames(dat_train), ]
prepf_test <- prepf[!rownames(prepf) %in% rownames(dat_train), ]

colnames(prepf_train) <- make.names(colnames(prepf_train))
colnames(prepf_test) <- make.names(colnames(prepf_test))

model <- train(DX ~ ., data = prepf_train, method = "glmnet", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$alpha == final_param$alpha & model$results$lambda == final_param$lambda, ]
model_final #good model overall

# try some tuning
tuning <- expand.grid(alpha = c(0.05, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 
                                0.12, 0.13, 0.14, 0.15, 0.2), 
                      lambda = c(0.01, 0.02, 0.03, 0.037, 0.038, 0.039, 0.04, 0.041, 0.042, 0.043,
                                 0.044, 0.045, 0.046, 0.047, 0.0475, 0.05, 0.06, 0.07, 0.08))
model2 <- train(DX ~ ., data = prepf_train, method = "glmnet", trControl = control, 
                tuneGrid = tuning)

final_param2 <- model2$bestTune #final model
model_final2 <- model2$results[model2$results$alpha == final_param2$alpha & model2$results$lambda == final_param2$lambda, ]
model_final2 #made it worse keep model 1

#look at top sites in model
top_vars <- varImp(model)$importance
top_vars$site <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$Overall, decreasing = T), ]
top_vars <- top_vars[top_vars$Overall != 0, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), Overall), y = Overall)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = Overall), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#look at results
predicted <- predict(model, newdata = prepf_test)
predicted

test_probs <- predict(model, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX)

#select training predictions using final hyperparamter values
select_idx <- model$pred$alpha == final_param$alpha & 
  model$pred$lambda == final_param$lambda
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_glmnet50diff <- model
prepf_glmnet50diff_testResults <- test_results

################ SVMLinear 50 features by differential testing ##################
model <- train(DX ~ ., data = prepf_train, method = "svmLinear", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$C == final_param$C, ]
model_final #good specificity and sensitivity

# try tuning
tuning <- expand.grid(C = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 0.95, 0.975, 
                            1, 10, 50, 100, 1000))
model2 <- train(DX ~ ., data = prepf_train, method = "svmLinear", trControl = control, 
                tuneGrid = tuning)

final_param2 <- model2$bestTune #final model
model_final2 <- model2$results[model2$results$C == final_param2$C, ]
model_final2 #worse; keep model 1

plot(model)

#look at results
predicted <- predict(model, newdata = prepf_test)

test_probs <- predict(model, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX)

#look at top sites in model
top_vars <- varImp(model)$importance
top_vars$site <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$PrePF, decreasing = T), ]
top_vars <- top_vars[top_vars$PrePF != 0, ] #994
ggplot(top_vars, aes(x = reorder(rownames(top_vars), PrePF), y = PrePF)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = PrePF), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#select training predictions using final hyperparamter values
select_idx <- model$pred$C == final_param$C
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_svml50diff <- model
prepf_svml50diff_testResults <- test_results

################ svmRadial 50 features by differential testing ###################
model <- train(DX ~ ., data = prepf_train, method = "svmRadial", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$sigma == final_param$sigma & model$results$C == final_param$C, ]
model_final #really good

# try tuning
tuning <- expand.grid(sigma = c(0.005, 0.0075, 0.008, 0.009, 0.0095, 0.01, 0.011, 0.012, 0.013,
                                0.014, 0.015, 0.016, 0.0175, 0.02, 0.05), 
                      C = c(0.01, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.75, 1, 2, 5, 
                            10, 15, 20, 25, 50))
model2 <- train(DX ~ ., data = prepf_train, method = "svmRadial", trControl = control, 
                tuneGrid = tuning)

final_param2 <- model2$bestTune #final model
model_final2 <- model2$results[model2$results$sigma == final_param2$sigma & model2$results$C == final_param2$C, ]
model_final2 #made it worse :(

plot(model)

#look at results
predicted <- predict(model, newdata = prepf_test)

test_probs <- predict(model, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX)

#look at top sites in model
top_vars <- varImp(model)$importance
top_vars$site <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$PrePF, decreasing = T), ]
top_vars <- top_vars[top_vars$PrePF != 0, ] #994
#top_vars <- top_vars[1:75, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), PrePF), y = PrePF)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = PrePF), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#select training predictions using final hyperparamter values
select_idx <- model$pred$C == final_param$C & model$pred$sigma == final_param$sigma
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_svmr50diff <- model
prepf_svmr50diff_testResults <- test_results

################ RF 50 features by differential testing ###############
model <- train(DX ~ ., data = prepf_train, method = "rf", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$mtry == final_param$mtry, ]
model_final #not very good

top_vars <- varImp(model)$importance
top_vars$gene <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$Overall, decreasing = T), ]
top_vars <- top_vars[top_vars$Overall != 0, ] #877
#top_vars <- top_vars[1:75, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), Overall), y = Overall)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = Overall), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

plot(model)

#look at results
predicted <- predict(model, newdata = prepf_test)
predicted

#also want predicted class probabilities 
test_probs <- predict(model, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX) 

#select training predictions using final hyperparamter values
select_idx <- model$pred$mtry == final_param$mtry
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_rf50diff <- model
prepf_rf50diff_testResults <- test_results

################ glmnet 15 features by differential testing ##################
sorted <- sort(col_var, decreasing = T, index.return = T)$ix[1:15]
sorted <- sorted + 1
dat_sub <- dat[, c(1,sorted)]
colnames(dat_sub) <- make.names(colnames(dat_sub))

#relevel for caret (takes first level as positive case)
dat_sub$DX <- as.factor(dat_sub$DX)
dat_sub$DX <- relevel(dat_sub$DX, ref = "PrePF")

#split data into training and test sets
#split 70:30 while maininting equivalent numbers of each DX in both groups
inTrain <- createDataPartition(dat_sub$DX, p = 0.7, list = F)

dat_train <- dat_sub[inTrain, ] #92
dat_test <- dat_sub[-inTrain, ] #38

prepf <- prepfvals[1:15, ] 

prepf <- dat[, colnames(dat) %in% rownames(prepf)] #all three
prepf <- cbind(DX = dat$DX, prepf)
#refactor for caret
prepf$DX <- as.factor(prepf$DX)
prepf$DX <- relevel(prepf$DX, ref = "PrePF")

prepf_train <- prepf[rownames(prepf) %in% rownames(dat_train), ]
prepf_test <- prepf[!rownames(prepf) %in% rownames(dat_train), ]

colnames(prepf_train) <- make.names(colnames(prepf_train))
colnames(prepf_test) <- make.names(colnames(prepf_test))

model <- train(DX ~ ., data = prepf_train, method = "glmnet", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$alpha == final_param$alpha & model$results$lambda == final_param$lambda, ]
model_final #surprisingly ok sensitivity for 15 features

# try some tuning
tuning <- expand.grid(alpha = c(0.01, 0.05, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11,
                                0.115, 0.12, 0.15, 0.2), 
                      lambda = c(0.001, 0.005, 0.0055, 0.006, 0.0065, 0.007, 0.00725, 0.0075, 
                                 0.00775, 0.008, 0.01))
model2 <- train(DX ~ ., data = prepf_train, method = "glmnet", trControl = control, 
                tuneGrid = tuning)

final_param2 <- model2$bestTune #final model
model_final2 <- model2$results[model2$results$alpha == final_param2$alpha & model2$results$lambda == final_param2$lambda, ]
model_final2 #better

plot(model2)

#look at top sites in model
top_vars <- varImp(model2)$importance
top_vars$site <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$Overall, decreasing = T), ]
top_vars <- top_vars[top_vars$Overall != 0, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), Overall), y = Overall)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = Overall), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#look at results
predicted <- predict(model2, newdata = prepf_test)
predicted

test_probs <- predict(model2, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX)

#select training predictions using final hyperparamter values
select_idx <- model2$pred$alpha == final_param2$alpha & 
  model2$pred$lambda == final_param2$lambda
g <- ggplot(model2$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model2$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_glmnet15diff <- model2
prepf_glmnet15diff_testResults <- test_results

################ SVMLinear 15 features by differential testing ##################
model <- train(DX ~ ., data = prepf_train, method = "svmLinear", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$C == final_param$C, ]
model_final #alright; okay sensitivity

# try tuning
tuning <- expand.grid(C = c(0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 10, 15, 20, 50, 100))
model2 <- train(DX ~ ., data = prepf_train, method = "svmLinear", trControl = control, 
                tuneGrid = tuning)

final_param2 <- model2$bestTune #final model
model_final2 <- model2$results[model2$results$C == final_param2$C, ]
model_final2 #improved all 3

plot(model2)

#look at results
predicted <- predict(model2, newdata = prepf_test)

test_probs <- predict(model2, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX)

#look at top sites in model
top_vars <- varImp(model2)$importance
top_vars$site <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$PrePF, decreasing = T), ]
top_vars <- top_vars[top_vars$PrePF != 0, ] #994
ggplot(top_vars, aes(x = reorder(rownames(top_vars), PrePF), y = PrePF)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = PrePF), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#select training predictions using final hyperparamter values
select_idx <- model2$pred$C == final_param2$C
g <- ggplot(model2$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model2$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_svml15diff <- model2
prepf_svml15diff_testResults <- test_results

################ svmRadial 15 features by differential testing ###################
model <- train(DX ~ ., data = prepf_train, method = "svmRadial", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$sigma == final_param$sigma & model$results$C == final_param$C, ]
model_final #pretty good sensitivity; but falling quickly

# try tuning
tuning <- expand.grid(sigma = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.055, 0.0575, 0.06, 0.0605, 0.061,
                                0.0615, 0.062, 0.0625, 0.063, 0.065, 0.07), 
                      C = c(0.01, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3, 4, 5, 10, 
                            15, 20, 25, 50, 100))
model2 <- train(DX ~ ., data = prepf_train, method = "svmRadial", trControl = control, 
                tuneGrid = tuning)

final_param2 <- model2$bestTune #final model
model_final2 <- model2$results[model2$results$sigma == final_param2$sigma & model2$results$C == final_param2$C, ]
model_final2 #improved sensitivity

plot(model2)
plot(model)

#look at results
predicted <- predict(model, newdata = prepf_test)
predicted2 <- predict(model2, newdata = prepf_test)

test_probs <- predict(model, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))
test_probs2 <- predict(model2, newdata = prepf_test, type = "prob")
test_results2 <- as.data.frame(cbind(predicted2, obs = prepf_test$DX, test_probs2))

confusionMatrix(data = predicted, reference = prepf_test$DX)
confusionMatrix(data = predicted2, reference = prepf_test$DX) #better so model 2

#look at top sites in model
top_vars <- varImp(model2)$importance
top_vars$site <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$PrePF, decreasing = T), ]
top_vars <- top_vars[top_vars$PrePF != 0, ] #994
#top_vars <- top_vars[1:75, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), PrePF), y = PrePF)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = PrePF), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#select training predictions using final hyperparamter values
select_idx <- model2$pred$C == final_param2$C & model2$pred$sigma == final_param2$sigma
g <- ggplot(model2$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model2$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

#plot ROC curves for our testing results
g <- ggplot(test_results2, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_svmr15diff <- model2
prepf_svmr15diff_testResults <- test_results2

################ RF 15 features by differential testing ###############
model <- train(DX ~ ., data = prepf_train, method = "rf", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$mtry == final_param$mtry, ]
model_final #getting better with less features but still bad

plot(model)

top_vars <- varImp(model)$importance
top_vars$gene <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$Overall, decreasing = T), ]
top_vars <- top_vars[top_vars$Overall != 0, ] #877
#top_vars <- top_vars[1:75, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), Overall), y = Overall)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = Overall), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#look at results
predicted <- predict(model, newdata = prepf_test)
predicted

#also want predicted class probabilities 
test_probs <- predict(model, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX) 

#select training predictions using final hyperparamter values
select_idx <- model$pred$mtry == final_param$mtry
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_rf15diff <- model
prepf_rf15diff_testResults <- test_results

################ glmnet 10 features by differential testing #################
sorted <- sort(col_var, decreasing = T, index.return = T)$ix[1:10]
sorted <- sorted + 1
dat_sub <- dat[, c(1,sorted)]
colnames(dat_sub) <- make.names(colnames(dat_sub))

#relevel for caret (takes first level as positive case)
dat_sub$DX <- as.factor(dat_sub$DX)
dat_sub$DX <- relevel(dat_sub$DX, ref = "PrePF")

#split data into training and test sets
#split 70:30 while maininting equivalent numbers of each DX in both groups
inTrain <- createDataPartition(dat_sub$DX, p = 0.7, list = F)

dat_train <- dat_sub[inTrain, ] #92
dat_test <- dat_sub[-inTrain, ] #38

prepf <- prepfvals[1:10, ] 

prepf <- dat[, colnames(dat) %in% rownames(prepf)] #all three
prepf <- cbind(DX = dat$DX, prepf)
#refactor for caret
prepf$DX <- as.factor(prepf$DX)
prepf$DX <- relevel(prepf$DX, ref = "PrePF")

prepf_train <- prepf[rownames(prepf) %in% rownames(dat_train), ]
prepf_test <- prepf[!rownames(prepf) %in% rownames(dat_train), ]

colnames(prepf_train) <- make.names(colnames(prepf_train))
colnames(prepf_test) <- make.names(colnames(prepf_test))

model <- train(DX ~ ., data = prepf_train, method = "glmnet", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$alpha == final_param$alpha & model$results$lambda == final_param$lambda, ]
model_final #ok sensitivity 

# try some tuning
tuning <- expand.grid(alpha = c(0.01, 0.05, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11,
                                0.115, 0.12, 0.15, 0.2), 
                      lambda = c(0.01, 0.02, 0.03, 0.04, 0.041, 0.042, 0.043, 0.044, 0.045, 0.046,
                                 0.05, 0.06))
model2 <- train(DX ~ ., data = prepf_train, method = "glmnet", trControl = control, 
                tuneGrid = tuning)

final_param2 <- model2$bestTune #final model
model_final2 <- model2$results[model2$results$alpha == final_param2$alpha & model2$results$lambda == final_param2$lambda, ]
model_final2 #better

plot(model2)

#look at top sites in model
top_vars <- varImp(model2)$importance
top_vars$site <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$Overall, decreasing = T), ]
top_vars <- top_vars[top_vars$Overall != 0, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), Overall), y = Overall)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = Overall), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#look at results
predicted <- predict(model2, newdata = prepf_test)
predicted

test_probs <- predict(model2, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX)

#select training predictions using final hyperparamter values
select_idx <- model2$pred$alpha == final_param2$alpha & 
  model2$pred$lambda == final_param2$lambda
g <- ggplot(model2$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model2$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_glmnet10diff <- model
prepf_glmnet10diff_testResults <- test_results

################ SVMLinear 10 features by differential testing ##################
model <- train(DX ~ ., data = prepf_train, method = "svmLinear", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$C == final_param$C, ]
model_final #alright; not the best sensitivity

# try tuning
tuning <- expand.grid(C = c(0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 10, 15, 20, 50, 100))
model2 <- train(DX ~ ., data = prepf_train, method = "svmLinear", trControl = control, 
                tuneGrid = tuning)

final_param2 <- model2$bestTune #final model
model_final2 <- model2$results[model2$results$C == final_param2$C, ]
model_final2 #big improvement in sensitivity; slight reduction in specificity

plot(model2)

#look at results
predicted <- predict(model2, newdata = prepf_test)

test_probs <- predict(model2, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX)

#look at top sites in model
top_vars <- varImp(model2)$importance
top_vars$site <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$PrePF, decreasing = T), ]
top_vars <- top_vars[top_vars$PrePF != 0, ] #994
ggplot(top_vars, aes(x = reorder(rownames(top_vars), PrePF), y = PrePF)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = PrePF), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#select training predictions using final hyperparamter values
select_idx <- model2$pred$C == final_param2$C
g <- ggplot(model2$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model2$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_svml10diff <- model2
prepf_svml10diff_testResults <- test_results

################ svmRadial 10 features by differential testing ###################
model <- train(DX ~ ., data = prepf_train, method = "svmRadial", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$sigma == final_param$sigma & model$results$C == final_param$C, ]
model_final #pretty decent for 10 features

# try tuning
tuning <- expand.grid(sigma = c(0.01, 0.025, 0.05, 0.06, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 
                                0.1), 
                      C = c(0.01, 0.1, 0.5, 1, 5, 10, 25, 50, 100, 1000))
model2 <- train(DX ~ ., data = prepf_train, method = "svmRadial", trControl = control, 
                tuneGrid = tuning)

final_param2 <- model2$bestTune #final model
model_final2 <- model2$results[model2$results$sigma == final_param2$sigma & model2$results$C == final_param2$C, ]
model_final2 #better

plot(model2)

#look at results
predicted <- predict(model2, newdata = prepf_test)

test_probs <- predict(model2, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX)

#look at top sites in model
top_vars <- varImp(model2)$importance
top_vars$site <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$PrePF, decreasing = T), ]
top_vars <- top_vars[top_vars$PrePF != 0, ] #994
#top_vars <- top_vars[1:75, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), PrePF), y = PrePF)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = PrePF), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

#select training predictions using final hyperparamter values
select_idx <- model2$pred$C == final_param2$C & model2$pred$sigma == final_param2$sigma
g <- ggplot(model2$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model2$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_svmr10diff <- model2
prepf_svmr10diff_testResults <- test_results

################ RF 10 features by differential testing ###############
model <- train(DX ~ ., data = prepf_train, method = "rf", trControl = control, tuneLength = 20)

final_param <- model$bestTune #final model
model_final <- model$results[model$results$mtry == final_param$mtry, ]
model_final #okay sensitivity

top_vars <- varImp(model)$importance
top_vars$gene <- rownames(top_vars)
top_vars <- top_vars[order(top_vars$Overall, decreasing = T), ]
top_vars <- top_vars[top_vars$Overall != 0, ] #877
#top_vars <- top_vars[1:75, ]
ggplot(top_vars, aes(x = reorder(rownames(top_vars), Overall), y = Overall)) + 
  geom_point(color = "blue", size = 4, alpha = 0.6) + 
  geom_segment(aes(x = rownames(top_vars), xend = rownames(top_vars), y = 0, yend = Overall), 
               color = "skyblue") + xlab("Variable") + ylab("Overall Importance") + 
  theme_light() + coord_flip()

plot(model)

#look at results
predicted <- predict(model, newdata = prepf_test)
predicted

#also want predicted class probabilities 
test_probs <- predict(model, newdata = prepf_test, type = "prob")
test_results <- as.data.frame(cbind(predicted, obs = prepf_test$DX, test_probs))

confusionMatrix(data = predicted, reference = prepf_test$DX) 

#select training predictions using final hyperparamter values
select_idx <- model$pred$mtry == final_param$mtry
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#make same plot with decision thresholds
g <- ggplot(model$pred[select_idx, ], aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) + geom_roc(n.cuts = 20) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))


#plot ROC curves for our testing results
g <- ggplot(test_results, aes(m = PrePF, d = factor(obs, levels = c("Cntrl", "PrePF")))) +
  geom_roc(n.cuts = 0) + coord_equal() + style_roc()
g + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round((calc_auc(g))$AUC, 4)))

prepf_rf10diff <- model
prepf_rf10diff_testResults <- test_results
################ xgboost 1000 features by differential testing ###############
library(xgboost)
library(readr)
library(stringr)
library(caret)
library(Ckmeans.1d.dp)
library(Matrix)
library(SHAPforxgboost)
library(DESeq2)
library(MLSeq)
library(dplyr)

sorted <- sort(col_var, decreasing = T, index.return = T)$ix[1:1000]
sorted <- sorted + 1
dat_sub <- dat[, c(1,sorted)]

#relevel for caret (takes first level as positive case)
dat_sub$DX <- as.factor(dat_sub$DX)
dat_sub$DX <- relevel(dat_sub$DX, ref = "PrePF")

#split data into training and test sets
#split 70:30 while maininting equivalent numbers of each DX in both groups
inTrain <- createDataPartition(dat_sub$DX, p = 0.7, list = F)

dat_train <- dat_sub[inTrain, ] #87
dat_test <- dat_sub[-inTrain, ] #36

prepf <- prepfvals[1:1000, ] 

prepf <- dat[, colnames(dat) %in% rownames(prepf)] #all three
prepf <- cbind(DX = dat$DX, prepf)
#refactor for caret
prepf$DX <- as.factor(prepf$DX)
prepf$DX <- relevel(prepf$DX, ref = "PrePF")

prepf_train <- prepf[rownames(prepf) %in% rownames(dat_train), ]
prepf_test <- prepf[!rownames(prepf) %in% rownames(dat_train), ]

colnames(prepf_train) <- make.names(colnames(prepf_train))
colnames(prepf_test) <- make.names(colnames(prepf_test))

train_labels <- (recode(prepf_train[ ,1], "PrePF" = 1, "Cntrl" = 0))
test_labels <- (recode(prepf_test[ ,1], "PrePF" = 1, "Cntrl" = 0))

#put data into xgb.DMatrix
dMat <- data.matrix(prepf_train[,2:ncol(prepf_train)])
dMat <- as(dMat, "dgCMatrix")
dMat2 <- data.matrix(prepf_test[,2:ncol(prepf_test)])
dMat2 <- as(dMat2, "dgCMatrix")

train <- list(data = dMat, label = train_labels)
test <- list(data = dMat2, label = test_labels)

dtrain <- xgb.DMatrix(data = train$data, label = train$label)
dtest <- xgb.DMatrix(data = test$data, label = test$label)

watchlist <- list(train = dtrain, test = dtest)
#best is depth 2, eta .3, nrounds 5
bst <- xgb.train(data = dtrain, max.depth = 4, eta = .4, nthread = 2, nrounds = 25,
                 watchlist = watchlist, objective = "binary:logistic",
                 eval.metric = "error")
pred <- predict(bst, test$data)
pred <- as.numeric(pred > .5)
pred
confusionMatrix(as.factor(pred), as.factor(test$label), positive = '1')

xgboost_1000diff <- bst
xgb.save(bst, "xgb_1000diff_bestModel")

importanceMatrix <- xgb.importance(feature_names = colnames(prepf_train)[-1],
                                   model = bst)
xgb.ggplot.importance(importanceMatrix)

#SHAP values
#need training data back into dataframe
train_df <- as.data.frame(as.matrix(dMat))
train_df <- as.matrix(train_df)
shap.plot.summary.wrap1(bst, X = train_df, top_n = 20)

#linear booster
bst3 <- xgb.train(data = dtrain, eta = .1, nthread = 2, nrounds = 20,
                  booster = "gblinear", lambda = 1, alpha = 0,
                  watchlist = watchlist, objective = "binary:logistic",
                  eval.metric = "error")
pred3 <- predict(bst3, test$data)
pred3 <- as.numeric(pred3 > .5)
pred3
confusionMatrix(as.factor(pred3), as.factor(test$label), positive = '1')
importanceMatrix <- xgb.importance(feature_names = colnames(prepf_train)[-1],
                                   model = bst3)
xgb.ggplot.importance(importanceMatrix[1:20,])
shap.plot.summary.wrap1(bst3, X = train_df, top_n = 20)

xgboost_1000diff_linear <- bst3
xgb.save(bst3, "xgb_1000diff_linear")

################# xgboost 100 features by differential testing ##################
sorted <- sort(col_var, decreasing = T, index.return = T)$ix[1:100]
sorted <- sorted + 1
dat_sub <- dat[, c(1,sorted)]

#relevel for caret (takes first level as positive case)
dat_sub$DX <- as.factor(dat_sub$DX)
dat_sub$DX <- relevel(dat_sub$DX, ref = "PrePF")

#split data into training and test sets
#split 70:30 while maininting equivalent numbers of each DX in both groups
inTrain <- createDataPartition(dat_sub$DX, p = 0.7, list = F)

dat_train <- dat_sub[inTrain, ] #87
dat_test <- dat_sub[-inTrain, ] #36

prepf <- prepfvals[1:100, ] 

prepf <- dat[, colnames(dat) %in% rownames(prepf)] #all three
prepf <- cbind(DX = dat$DX, prepf)
#refactor for caret
prepf$DX <- as.factor(prepf$DX)
prepf$DX <- relevel(prepf$DX, ref = "PrePF")

prepf_train <- prepf[rownames(prepf) %in% rownames(dat_train), ]
prepf_test <- prepf[!rownames(prepf) %in% rownames(dat_train), ]

colnames(prepf_train) <- make.names(colnames(prepf_train))
colnames(prepf_test) <- make.names(colnames(prepf_test))

train_labels <- (recode(prepf_train[ ,1], "PrePF" = 1, "Cntrl" = 0))
test_labels <- (recode(prepf_test[ ,1], "PrePF" = 1, "Cntrl" = 0))

#put data into xgb.DMatrix
dMat <- data.matrix(prepf_train[,2:ncol(prepf_train)])
dMat <- as(dMat, "dgCMatrix")
dMat2 <- data.matrix(prepf_test[,2:ncol(prepf_test)])
dMat2 <- as(dMat2, "dgCMatrix")

train <- list(data = dMat, label = train_labels)
test <- list(data = dMat2, label = test_labels)

dtrain <- xgb.DMatrix(data = train$data, label = train$label)
dtest <- xgb.DMatrix(data = test$data, label = test$label)

watchlist <- list(train = dtrain, test = dtest)
#best is depth 2, eta .1, 20 rounds
bst <- xgb.train(data = dtrain, max.depth = 2, eta = .1, nthread = 2, nrounds = 20,
                 watchlist = watchlist, objective = "binary:logistic",
                 eval.metric = "error")
pred <- predict(bst, test$data)
pred <- as.numeric(pred > .5)
pred
confusionMatrix(as.factor(pred), as.factor(test$label), positive = '1')

xgboost_100diff <- bst
xgb.save(bst, "xgb_100diff_bestModel")

importanceMatrix <- xgb.importance(feature_names = colnames(prepf_train)[-1],
                                   model = bst)
xgb.ggplot.importance(importanceMatrix)

#SHAP values
#need training data back into dataframe
train_df <- as.data.frame(as.matrix(dMat))
train_df <- as.matrix(train_df)
shap.plot.summary.wrap1(bst, X = train_df, top_n = 14)

#linear booster
bst3 <- xgb.train(data = dtrain, eta = .2, nthread = 2, nrounds = 100,
                  booster = "gblinear", lambda = 1, alpha = 0,
                  watchlist = watchlist, objective = "binary:logistic",
                  eval.metric = "error")
pred3 <- predict(bst3, test$data)
pred3 <- as.numeric(pred3 > .5)
pred3
confusionMatrix(as.factor(pred3), as.factor(test$label), positive = '1')
importanceMatrix <- xgb.importance(feature_names = colnames(prepf_train)[-1],
                                   model = bst3)
xgb.ggplot.importance(importanceMatrix[1:20,])
shap.plot.summary.wrap1(bst3, X = train_df, top_n = 20)

xgboost_linear_100diff <- bst3

#################### xgboost 50 features by differential testing ###############
sorted <- sort(col_var, decreasing = T, index.return = T)$ix[1:50]
sorted <- sorted + 1
dat_sub <- dat[, c(1,sorted)]

#relevel for caret (takes first level as positive case)
dat_sub$DX <- as.factor(dat_sub$DX)
dat_sub$DX <- relevel(dat_sub$DX, ref = "PrePF")

#split data into training and test sets
#split 70:30 while maininting equivalent numbers of each DX in both groups
inTrain <- createDataPartition(dat_sub$DX, p = 0.7, list = F)

dat_train <- dat_sub[inTrain, ] #87
dat_test <- dat_sub[-inTrain, ] #36

prepf <- prepfvals[1:50, ] 

prepf <- dat[, colnames(dat) %in% rownames(prepf)] #all three
prepf <- cbind(DX = dat$DX, prepf)
#refactor for caret
prepf$DX <- as.factor(prepf$DX)
prepf$DX <- relevel(prepf$DX, ref = "PrePF")

prepf_train <- prepf[rownames(prepf) %in% rownames(dat_train), ]
prepf_test <- prepf[!rownames(prepf) %in% rownames(dat_train), ]

colnames(prepf_train) <- make.names(colnames(prepf_train))
colnames(prepf_test) <- make.names(colnames(prepf_test))

train_labels <- (recode(prepf_train[ ,1], "PrePF" = 1, "Cntrl" = 0))
test_labels <- (recode(prepf_test[ ,1], "PrePF" = 1, "Cntrl" = 0))

#put data into xgb.DMatrix
dMat <- data.matrix(prepf_train[,2:ncol(prepf_train)])
dMat <- as(dMat, "dgCMatrix")
dMat2 <- data.matrix(prepf_test[,2:ncol(prepf_test)])
dMat2 <- as(dMat2, "dgCMatrix")

train <- list(data = dMat, label = train_labels)
test <- list(data = dMat2, label = test_labels)

dtrain <- xgb.DMatrix(data = train$data, label = train$label)
dtest <- xgb.DMatrix(data = test$data, label = test$label)

watchlist <- list(train = dtrain, test = dtest)
#best is depth 4, eta .4, nrounds 10
bst <- xgb.train(data = dtrain, max.depth = 4, eta = .4, nthread = 2, nrounds = 10,
                 watchlist = watchlist, objective = "binary:logistic",
                 eval.metric = "error")
pred <- predict(bst, test$data)
pred <- as.numeric(pred > .5)
pred
confusionMatrix(as.factor(pred), as.factor(test$label), positive = '1')

xgboost_50diff <- bst
xgb.save(bst, "xgb_50diff_bestModel")

importanceMatrix <- xgb.importance(feature_names = colnames(prepf_train)[-1],
                                   model = bst)
xgb.ggplot.importance(importanceMatrix)

#SHAP values
#need training data back into dataframe
train_df <- as.data.frame(as.matrix(dMat))
train_df <- as.matrix(train_df)
shap.plot.summary.wrap1(bst, X = train_df, top_n = 20)

#linear booster
bst3 <- xgb.train(data = dtrain, eta = .1, nthread = 2, nrounds = 25,
                  booster = "gblinear", lambda = 1, alpha = 0,
                  watchlist = watchlist, objective = "binary:logistic",
                  eval.metric = "error")
pred3 <- predict(bst3, test$data)
pred3 <- as.numeric(pred3 > .5)
pred3
confusionMatrix(as.factor(pred3), as.factor(test$label), positive = '1')
importanceMatrix <- xgb.importance(feature_names = colnames(prepf_train)[-1],
                                   model = bst3)
xgb.ggplot.importance(importanceMatrix[1:20,])
shap.plot.summary.wrap1(bst3, X = train_df, top_n = 20)

##################### xgboost 15 features by differential testing ############
sorted <- sort(col_var, decreasing = T, index.return = T)$ix[1:15]
sorted <- sorted + 1
dat_sub <- dat[, c(1,sorted)]

#relevel for caret (takes first level as positive case)
dat_sub$DX <- as.factor(dat_sub$DX)
dat_sub$DX <- relevel(dat_sub$DX, ref = "PrePF")

#split data into training and test sets
#split 70:30 while maininting equivalent numbers of each DX in both groups
inTrain <- createDataPartition(dat_sub$DX, p = 0.7, list = F)

dat_train <- dat_sub[inTrain, ] #87
dat_test <- dat_sub[-inTrain, ] #36

prepf <- prepfvals[1:15, ] 

prepf <- dat[, colnames(dat) %in% rownames(prepf)] #all three
prepf <- cbind(DX = dat$DX, prepf)
#refactor for caret
prepf$DX <- as.factor(prepf$DX)
prepf$DX <- relevel(prepf$DX, ref = "PrePF")

prepf_train <- prepf[rownames(prepf) %in% rownames(dat_train), ]
prepf_test <- prepf[!rownames(prepf) %in% rownames(dat_train), ]

colnames(prepf_train) <- make.names(colnames(prepf_train))
colnames(prepf_test) <- make.names(colnames(prepf_test))

train_labels <- (recode(prepf_train[ ,1], "PrePF" = 1, "Cntrl" = 0))
test_labels <- (recode(prepf_test[ ,1], "PrePF" = 1, "Cntrl" = 0))

#put data into xgb.DMatrix
dMat <- data.matrix(prepf_train[,2:ncol(prepf_train)])
dMat <- as(dMat, "dgCMatrix")
dMat2 <- data.matrix(prepf_test[,2:ncol(prepf_test)])
dMat2 <- as(dMat2, "dgCMatrix")

train <- list(data = dMat, label = train_labels)
test <- list(data = dMat2, label = test_labels)

dtrain <- xgb.DMatrix(data = train$data, label = train$label)
dtest <- xgb.DMatrix(data = test$data, label = test$label)

watchlist <- list(train = dtrain, test = dtest)
#best is depth 3, eta .4, nrounds 15
bst <- xgb.train(data = dtrain, max.depth = 3, eta = .4, nthread = 2, nrounds = 15,
                 watchlist = watchlist, objective = "binary:logistic",
                 eval.metric = "error")
pred <- predict(bst, test$data)
pred <- as.numeric(pred > .5)
pred
confusionMatrix(as.factor(pred), as.factor(test$label), positive = '1')

xgboost_15diff <- bst
xgb.save(bst, "xgb_15diff_bestModel")

importanceMatrix <- xgb.importance(feature_names = colnames(prepf_train)[-1],
                                   model = bst)
xgb.ggplot.importance(importanceMatrix)

#SHAP values
#need training data back into dataframe
train_df <- as.data.frame(as.matrix(dMat))
train_df <- as.matrix(train_df)
shap.plot.summary.wrap1(bst, X = train_df)

################### feature comparison #####################
glmnet_100diff_pred <- predictors(glmnet_100diff) #72
glmnet_50diff_pred <- predictors(glmnet_50diff) #30
rf_50diff_pred <- predictors(rf_50diff) #50
svmr_100diff_pred <- predictors(svmr_100diff) #100
svmr_50diff_pred <- predictors(svmr_50diff) #50
svmr_15diff_pred <- predictors(svmr_15diff) #15
xgb_1000var_pred <- head(prepf_imp, 50)$Feature

dataList <- list(glmnet_100Diff = glmnet_100diff_pred,
                 glmnet_50Diff = glmnet_50diff_pred,
                 rf_50Diff = rf_50diff_pred,
                 svmr_100Diff = svmr_100diff_pred,
                 svmr_50Diff = svmr_50diff_pred,
                 svmr_15Diff = svmr_15diff_pred,
                 xgb_1000Var = xgb_1000var_pred)
dataList100s <- list(glmnet_100Diff = glmnet_100diff_pred,
                     svmr_100Diff = svmr_100diff_pred)
dataList50s <- list(glmnet_50Diff = glmnet_50diff_pred,
                    rf_50Diff = rf_50diff_pred,
                    svmr_50Diff = svmr_50diff_pred,
                    xgb_1000Var = xgb_1000var_pred)

v = Venn(dataList)
v.100 <- Venn(dataList100s)
v.50 <- Venn(dataList50s)
setmap(v, element_clustering = F, set_clustering = F, element_fontsize = 4)
setmap(v.100, element_clustering = F, set_clustering = F, element_fontsize = 6)
setmap(v.50, element_clustering = F, set_clustering = F, element_fontsize = 6)
