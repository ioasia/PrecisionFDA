### PrecisionFDA Challenge

## Libraries

## Libraries

## Libraries

# install.packages('easypackages')
library(easypackages)

toInstall <- c("reshape2", "ggplot2", "circlize", 'mixOmics', 'ROCR', 'ggrepel')
# install.packages(toInstall)

# source("https://bioconductor.org/biocLite.R")
toInstallBio <- c('biomaRt', 'pcaMethods')
# biocLite(toInstallBio)

libraries(c(toInstall, toInstallBio))

# Data
prot <- read.table("~/Desktop/PrecisionFDA/Files/MissForestImputedData.txt")
clin <- read.table('~/Desktop/PrecisionFDA/Files/train_cli.tsv', skip = 1)
clas <- read.table('~/Desktop/PrecisionFDA/Files/sum_tab_1.csv', skip = 1, sep = ',')

# Log2 transform
prot <- as.matrix(log2(prot))

# Remove mislabeled samples
colsToKeep <- clas$V2 !=1

# Find features for predicting MSI
Y <- as.vector(t(clin$V3[colsToKeep]))
# table(Y)

X <- t(prot[, colsToKeep])

# Sparse partial least squares for MSI 
splsdaRes <- splsda(X, Y, ncomp = 10)

# Evaluate model with leave-one-out
perfPlsda <- perf(splsdaRes, validation = 'loo', 
                  progressBar = FALSE, auc = TRUE) 

# Error rates
plot(perfPlsda, 
     col = color.mixo(1:3), 
     sd = TRUE, legend.position = "horizontal")

# Tuning sPLS-DA
# Grid of possible keepX values that will be tested for each comp 
list.keepX <- seq(10, 200, 20)

tunePlsda <- tune.splsda(X, Y, ncomp = 2, validation = 'loo', 
                         progressBar = TRUE, dist = 'max.dist',
                         test.keepX = list.keepX)


# tune.splsda.srbct  #the various outputs
choice.keepX <- tunePlsda$choice.keepX[1:2]

## sPLS-DA function
splsdaRes <- splsda(X, Y, ncomp = 2, keepX = choice.keepX)

## Plot components
p <- plotIndiv(splsdaRes, ind.names = Y, comp = c(1,2), ellipse = TRUE, legend = TRUE, title = 'SPLS-DA Components', style = 'ggplot2', legend.title = NULL)

## Keep proteins of 1st component
selectedFeatures <- rownames(selectVar(splsdaRes, comp = 1)$value)
selectedFeatures.2 <- rownames(selectVar(splsdaRes, comp = 2)$value)
selectedFeatures <- as.character(unique(c(selectedFeatures,selectedFeatures.2)))

library(caret)
library(e1071)
# get selected features from prot df 
selectedFeature.df <- prot[selectedFeatures,]

# get samples with the correctly classified
sel.fea.col.df <- data.frame(t(selectedFeature.df[, colsToKeep]))

# add the msi class information
sel.fea.col.df$class <- factor(t(clin$V3[colsToKeep]))

#prepare test data which are the all misclassified samples
selected.mis.cls <- data.frame(t(selectedFeature.df[, !colsToKeep]))
selected.mis.cls$class <- factor(t(clin$V3[ ! colsToKeep]))
testing.df <- selected.mis.cls

#use all 68 sample for classification
training.df <- sel.fea.col.df
levels(training.df$class) <- make.names(levels(factor(training.df$class)))
levels(testing.df$class) <- make.names(levels(factor(testing.df$class)))

ctrl <- trainControl(method = "repeatedcv", number = 5, savePredictions = TRUE, classProbs =  TRUE)

SVMgrid <- expand.grid(sigma = c(0.0010,0.0015,0.002,0.022,0.025,0.027, 0.03, 0.035), C = c(0.5,1,1.5,1.65,2,2.2,2.44,2.6,2.8,3))
set.seed(10)
mod.fit.msi <- train(class~., data=training.df, method="svmRadial", trControl = ctrl,
                 tuneGrid = SVMgrid, preProc = c("scale","YeoJohnson"), verbose=FALSE)

test.label <- read.table("~/Desktop/PrecisionFDA/Files/test_cli.txt", header = T)
levels(test.label$msi) <- make.names(levels(factor(test.label$msi)))

test.df <- read.table("~/Desktop/PrecisionFDA/Files/MissForestImputedTESTData.txt", header = T)
test.df <- t(log2(as.matrix(test.df[selectedFeatures,])))

pred.msi <- predict(mod.fit.msi, test.df)
pred.msi.prob <- predict(mod.fit.msi, test.df, type = "prob")
confusionMatrix(data=pred.msi, factor(test.label$msi))

pred.msi.df <- data.frame(Observation.Msi = test.label$msi, Prediction.msi = pred.msi, pred.msi.prob)
rownames(pred.msi.df) <- rownames(test.label)

msi.pred <- write.table(pred.msi.df,"~/Desktop/PrecisionFDA/msi.test.pred.txt", row.names = T, sep = "\t", quote = FALSE)

################ building classifier for geneder
library(easypackages)

toInstall <- c("reshape2", "ggplot2", "circlize", 'mixOmics', 'ROCR', 'ggrepel')
# install.packages(toInstall)

# source("https://bioconductor.org/biocLite.R")
toInstallBio <- c('biomaRt', 'pcaMethods')
# biocLite(toInstallBio)

libraries(c(toInstall, toInstallBio))

# Data
prot <- read.table("~/Desktop/PrecisionFDA/Files/MissForestImputedData.txt")
clin <- read.table('~/Desktop/PrecisionFDA/Files/train_cli.tsv', skip = 1)
clas <- read.table('~/Desktop/PrecisionFDA/Files/sum_tab_1.csv', skip = 1, sep = ',')

# Log2 transform
prot <- as.matrix(log2(prot))

# Remove mislabeled samples
colsToKeep <- clas$V2 !=1

# Find features for predicting MSI
Y <- as.vector(t(clin$V2[colsToKeep]))
# table(Y)

X <- t(prot[, colsToKeep])

# Sparse partial least squares for MSI 
splsdaRes <- splsda(X, Y, ncomp = 10)

# Evaluate model with leave-one-out
perfPlsda <- perf(splsdaRes, validation = 'loo', 
                  progressBar = FALSE, auc = TRUE) 

# Error rates
plot(perfPlsda, 
     col = color.mixo(1:3), 
     sd = TRUE, legend.position = "horizontal")

# Tuning sPLS-DA
# Grid of possible keepX values that will be tested for each comp 
list.keepX <- seq(10, 200, 10)

tunePlsda <- tune.splsda(X, Y, ncomp = 2, validation = 'loo', 
                         progressBar = TRUE, dist = 'max.dist',
                         test.keepX = list.keepX)


# tune.splsda.srbct  #the various outputs
choice.keepX <- tunePlsda$choice.keepX[1:2]

## sPLS-DA function
splsdaRes <- splsda(X, Y, ncomp = 2, keepX = choice.keepX)

## Plot components
p <- plotIndiv(splsdaRes, ind.names = Y, comp = c(1,2), ellipse = TRUE, legend = TRUE, title = 'SPLS-DA Components', style = 'ggplot2', legend.title = NULL)

## Keep proteins of 1st component
selectedFeatures <- rownames(selectVar(splsdaRes, comp = 1)$value)
selectedFeatures.2 <- rownames(selectVar(splsdaRes, comp = 2)$value)
#selectedFeatures.3 <- rownames(selectVar(splsdaRes, comp = 3)$value)
#selectedFeatures.4 <- rownames(selectVar(splsdaRes, comp = 4)$value)
selectedFeatures <- as.character(unique(c(selectedFeatures,selectedFeatures.2)))

library(caret)
library(e1071)
# get selected features from prot df 
selectedFeature.df <- prot[selectedFeatures,]

# get samples with the correctly classified
sel.fea.col.df <- data.frame(t(selectedFeature.df[, colsToKeep]))

# add the msi class information
sel.fea.col.df$class <- factor(t(clin$V2[colsToKeep]))

#prepare test data which are the all misclassified samples
selected.mis.cls <- data.frame(t(selectedFeature.df[, !colsToKeep]))
selected.mis.cls$class <- factor(t(clin$V2[ ! colsToKeep]))
testing.df <- selected.mis.cls

#use all 68 sample for classification
training.df <- sel.fea.col.df
levels(training.df$class) <- make.names(levels(factor(training.df$class)))
levels(testing.df$class) <- make.names(levels(factor(testing.df$class)))

ctrl <- trainControl(method = "repeatedcv", number = 5, savePredictions = TRUE, classProbs =  TRUE)

SVMgrid <- expand.grid(sigma = c(0.0010,0.0015,0.002,0.022,0.025,0.027, 0.03, 0.035), C = c(0.5,1,1.5,1.65,2,2.2,2.44,2.6,2.8,3))
set.seed(10)
mod.fit.gender <- train(class~., data=training.df, method="svmRadial", trControl = ctrl,
                     tuneGrid = SVMgrid, preProc = c("scale","YeoJohnson"), verbose=FALSE)

pred = predict(mod.fit.gender, newdata=testing.df)
pred.prob <-predict(mod.fit.gender, newdata = testing.df, type = "prob")
confusionMatrix(data=pred, factor(testing.df$class))


test.label <- read.table("~/Desktop/PrecisionFDA/Files/test_cli.txt", header = T)
levels(test.label$gender) <- make.names(levels(factor(test.label$gender)))

test.df <- read.table("~/Desktop/PrecisionFDA/Files/MissForestImputedTESTData.txt", header = T)
test.df <- t(log2(as.matrix(test.df[selectedFeatures,])))

pred.gender <- predict(mod.fit.gender, test.df)
pred.gender.prob <- predict(mod.fit.gender, test.df, type = "prob")
confusionMatrix(data=pred.gender, factor(test.label$gender))

pred.gender.df <- data.frame(Observation.Gender = test.label$gender, Prediction.Gender = pred.gender, pred.gender.prob)
rownames(pred.gender.df) <- test.label$sample

pred.gender.df$max.prob <- apply(pred.gender.df, 1, function(x) max(as.numeric(x[3]), as.numeric(x[4])))


gender.pred <- write.table(pred.gender.df,"~/Desktop/PrecisionFDA/gender.test.pred.txt", row.names = T, sep = "\t", quote = FALSE)

##### merge two results
msi.df <- read.table("~/Desktop/PrecisionFDA/msi.test.pred.txt", header = T)
msi.df$MSI.PRED <- ifelse(msi.df$Observation.Msi == msi.df$Prediction.msi, 0, 1)
gender.df <- read.table("~/Desktop/PrecisionFDA/gender.test.pred.txt", header = T)
gender.df$GENDER.PRED <- ifelse(gender.df$Observation.Gender == gender.df$Prediction.Gender, 0, 1)
gender.df[gender.df$max.prob < 0.9,][6] <- as.numeric(0)
gender.df <- gender.df[,-5]

final.pred <- data.frame(msi.df, gender.df)
final.pred$Merged.Pred <- as.numeric(final.pred$MSI.PRED) + as.numeric(final.pred$GENDER.PRED)
final.pred$Merged.Pred <- ifelse(final.pred$Merged.Pred > 0, 1, 0)
rownames(final.pred) <- test.label$sample

df <- data.frame(sample = test.label$sample, mismatch = final.pred$Merged.Pred)
write.csv(df, "~/Desktop/PrecisionFDA/Files/subchallenge_1.csv",quote = FALSE, sep = ",", row.names = FALSE)

