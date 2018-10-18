### PrecisionFDA Challenge

## Libraries

# install.packages('easypackages')
library(easypackages)

toInstall <- c("reshape2", "ggplot2", "circlize", 'mixOmics', 'ROCR', 'ggrepel', 'DMwR')
# install.packages(toInstall)

# source("https://bioconductor.org/biocLite.R")
toInstallBio <- c('biomaRt', 'pcaMethods')
# biocLite(toInstallBio)

libraries(c(toInstall, toInstallBio))


# Data
prot <- read.table('~/Desktop/PrecisionFDAChallenge/train_pro.tsv')
clin <- read.table('~/Desktop/PrecisionFDAChallenge/train_cli.tsv', skip = 1)
clas <- read.table('~/Desktop/PrecisionFDAChallenge/sum_tab_1.csv', skip = 1, sep = ',')


## Running parameters
transf <- 'Unlogged'
imput <- 'RandomForest'
featRed <- 'SPLSDA'
classif <- 'SVM'
augm <- 'None'
comp <- 'Gender_MSI'

# set.seed(123)

# Log2 transform
# prot <- log2(prot)


# Impute via local least squares imputation (Valikangas et. al. paper)

# Remove missing values in 30% of samples
# toKeep <- apply(prot, 1, function(i) sum(is.na(i)) < round(ncol(prot)/3))
# prot <- prot[toKeep, ]
# 
# # Impute
# protImputed <- llsImpute(t(prot), k = 150, correlation="spearman", allVariables=TRUE)
# prot <- t(completeObs(protImputed))

# Impute via random forest
prot <- read.table('imputation_result/MissForestImputedData.txt')


# Remove mislabeled samples
colsToKeep <- clas$V2 !=1


# Find features for predicting MSI
# Y <- as.vector(t(clin$V3[colsToKeep]))


# Include 4 groups (MSI and Gender)
Y <- as.factor(paste(clin$V2[colsToKeep], clin$V3[colsToKeep]))

table(Y)

X <- t(prot[, colsToKeep])


# Train Gender in MSI separately
# Y_all <- as.factor(paste(clin$V2[colsToKeep], clin$V3[colsToKeep]))
# Y_gender_high <- Y_all %in% c('Female MSI-High', 'Male MSI-High')
# Y_gender_low <- Y_all %in% c('Female MSI-Low/MSS', 'Male MSI-Low/MSS')

# X <- t(prot[, colsToKeep][,Y_gender_high])
# Y <- as.factor(paste(clin$V2[colsToKeep][Y_gender_high], 
#                      clin$V3[colsToKeep][Y_gender_high]))


## Data augmentation
# smoteInput <- cbind.data.frame(X,Y)
# smoteX <- SMOTE(Y~., smoteInput,  perc.over = 500, perc.under= 500)
# table(smoteX$Y)

## Remove introduced NAs 
# idxNa <- complete.cases(smoteX)
# smoteX[!idxNa, ncol(smoteX)]
# 
# Y <- smoteX[idxNa, ncol(smoteX)]
# X <- smoteX[idxNa, -ncol(smoteX)]



# Sparse partial least squares for MSI 
splsdaRes <- splsda(X, Y, ncomp = 10)

# Evaluate model with leave-one-out
perfPlsda <- perf(splsdaRes, validation = 'loo',
progressBar = FALSE, auc = TRUE)

# Evaluate model with cross - validation
# perfPlsda <- perf(splsdaRes, validation =  'Mfold', folds = 10,
#                   progressBar = FALSE, auc = TRUE)


# Error rates -This plot shows that the error rate reaches minimum at 2 Components
pdf(paste0('Feature_selection_classification_error', '_', paste(transf, imput, featRed, classif, augm, comp, sep = '_'),  '.pdf'), width = 8, height = 8)
plot(perfPlsda, 
     col = color.mixo(1:3), 
     sd = TRUE, legend.position = "horizontal")
dev.off()


# Tuning sPLS-DA - Choice of components based on classification error plot
nComps <- 2
# Grid of possible keepX values that will be tested for each component (This step tests for different number of proteins to include in its component)
list.keepX <- seq(10, 200, 20)

tunePlsda <- tune.splsda(X, Y, ncomp = nComps, validation = 'loo',
progressBar = TRUE, dist = 'max.dist',
test.keepX = list.keepX)

# tunePlsda <- tune.splsda(X, Y, ncomp = nComps, validation = 'Mfold', folds = 10,
#                          progressBar = TRUE, dist = 'max.dist',
#                          test.keepX = list.keepX)



# tune.splsda.srbct  #the various outputs (10 for component 1, 130 for component 2)

choice.keepX <- tunePlsda$choice.keepX[1:nComps]

## sPLS-DA function
splsdaRes <- splsda(X, Y, ncomp = 2, keepX = choice.keepX)


## Plot components
comp1 = 1
comp2 = 2
p <- plotIndiv(splsdaRes, ind.names = Y, comp = c(comp1, comp2), ellipse = TRUE, legend = TRUE,legend.position = 'bottom', title = 'SPLS-DA Components', style = 'ggplot2', legend.title = NULL)

pdf(paste0('BC_SPLS_Components_',comp1,'_', comp2, '_', paste(transf, imput, featRed, classif, augm, comp, sep = '_'), '.pdf'), width = 8, height = 8)
plot(p$graph)
dev.off()


## Keep proteins of components
selectedFeatures <- unlist(sapply(1:nComps, function(i) rownames(selectVar(splsdaRes, comp = i)$value)))



## SVM classifier based on first two principal components
library(caret)

# prComp <- prcomp(t(prot[selectedFeatures, ]))

data_part <- createDataPartition(y = paste(clin$V2[colsToKeep], clin$V3[colsToKeep]),
                                 p = 0.7, list = F)

# genderHigh <- grepl('MSI-High', Y_all)
# Y_all_gender_high <- Y_all[genderHigh]

# X_all_gender_high <-  t(prot[, colsToKeep][,Y_gender_high])

# data_part <- createDataPartition(Y_all_gender_high,
#                                  p = 0.7, list = F)

training <- data.frame(feat = t(prot[selectedFeatures, colsToKeep][, data_part]),
                       class =  paste(clin$V2, clin$V3)[colsToKeep][data_part])

# training <- data.frame(feat = X_all_gender_high[data_part, selectedFeatures],
# class =  Y_all_gender_high[data_part])
# training <- training[complete.cases(training), ]

testing <- data.frame(feat =  t(prot[selectedFeatures, colsToKeep][, -data_part]),
                      class =  paste(clin$V2, clin$V3)[colsToKeep][-data_part])

# testing <- data.frame(feat = X_all_gender_high[-data_part,selectedFeatures],
#                        class =  Y_all_gender_high[-data_part])

ctrl <- trainControl(method = "repeatedcv", number = 5, savePredictions = TRUE)

SVMgrid <- expand.grid(sigma = c(0.05,0.0456,0.0577), C = c(1.5,1.596,1.65,1.89,1.95,2,2.2,2.44))

mod_fit <- train(class~., data=training, method="svmRadial", trControl = ctrl,
                     tuneGrid = SVMgrid, preProc = c("scale","YeoJohnson"), verbose=FALSE)

# print(mod_fit)
pred = predict(mod_fit, newdata=testing)
confusionMatrix(data=pred, testing$class)



## Predict mislabelled 
misLabelled <-  data.frame(feat = t(prot[selectedFeatures, !colsToKeep]),
                           class =  paste(clin$V2, clin$V3)[!colsToKeep])

sum(predict(mod_fit, newdata = misLabelled) !=  paste(clin$V2, clin$V3)[!colsToKeep])/ nrow(misLabelled)

idx <- which(predict(mod_fit, newdata = misLabelled) == paste(clin$V2, clin$V3)[!colsToKeep])



## Vizualize classifier
# library(e1071)
# 
# model = svm(class~., data=training)
# plot(model, testing)


# PCA for vizualization
prComp <- prcomp(t(prot[selectedFeatures, ]))

## First two views
explainedVar <- round(prComp$sdev^2/sum(prComp$sdev^2), 2) * 100

## Plot
data <- data.frame(obsnames=row.names(prComp$x), prComp$x, Type = clin$V3)

pdf(paste0('Principal_Component_misclassified_', paste(transf, imput, featRed, classif,augm, comp, sep = '_'), '.pdf'), width = 6, height = 6)
p <- ggplot(data, aes(x=PC1, y=PC2)) + 
  geom_point(size = 5, aes(colour = Type)) + 
  geom_point(data =  data[!colsToKeep, ],
             shape = 21, size = 8, colour = 'red', stroke = 2) + 
  geom_label_repel(data = data[as.character(clin$V1[!colsToKeep][idx]), ], 
                   label = 'missed', nudge_y = 0.5) + 
  geom_text(aes(label = sapply(strsplit(rownames(data), "_"), function(i) i[[2]]))) + 
  scale_x_continuous(name = paste0("PC1 ","(", explainedVar[1], '%)')) +
  scale_y_continuous(name = paste0("PC2 ","(", explainedVar[2], '%)')) + 
  theme_bw() + 
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5), 
        text = element_text(size = 10), legend.text = element_text(size = 5), 
        legend.title = element_text(size = 8)) + 
  guides(colour = guide_legend(title = 'Subtype', nrow=1)) 

p + coord_flip()
dev.off()

###