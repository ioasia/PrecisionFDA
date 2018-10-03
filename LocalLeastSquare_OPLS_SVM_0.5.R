### PrecisionFDA Challenge

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
prot <- read.table('~/Desktop/PrecisionFDAChallenge/train_pro.tsv')
clin <- read.table('~/Desktop/PrecisionFDAChallenge/train_cli.tsv', skip = 1)
clas <- read.table('~/Desktop/PrecisionFDAChallenge/sum_tab_1.csv', skip = 1, sep = ',')



# Log2 transform
prot <- log2(prot)


# Impute via local least squares imputation (Valikangas et. al. paper)

# Remove missing values in 30% of samples 
toKeep <- apply(prot, 1, function(i) sum(is.na(i)) < round(ncol(prot)/3))
prot <- prot[toKeep, ]

# Impute
protImputed <- llsImpute(t(prot), k = 150, correlation="spearman", allVariables=TRUE)
prot <- t(completeObs(protImputed))


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

# pdf(paste0('BC_SPLS_Components_',comp1,'_', comp2, '.pdf'), width = 5, height = 5)
plot(p$graph)
# dev.off()


## Keep proteins of 1st component
selectedFeatures <- rownames(selectVar(splsdaRes, comp = 1)$value)



## SVM classifier based on first two principal components
library(caret)

prComp <- prcomp(t(prot[selectedFeatures, ]))

data_part <- createDataPartition(y = clin$V3[colsToKeep], 
                                 p = 0.7, list = F)

training <- data.frame(PC1 = prComp$x[,1][colsToKeep][data_part], 
                       PC2 = prComp$x[,2][colsToKeep][data_part],
                       class =  clin$V3[colsToKeep][data_part])

testing <- data.frame(PC1 =  prComp$x[,1][colsToKeep][-data_part], 
                      PC2 =  prComp$x[,2][colsToKeep][-data_part],
                      class =  clin$V3[colsToKeep][-data_part])


ctrl <- trainControl(method = "repeatedcv", number = 5, savePredictions = TRUE)

SVMgrid <- expand.grid(sigma = c(0.05,0.0456,0.0577), C = c(1.5,1.596,1.65,1.89,1.95,2,2.2,2.44))

mod_fit <- train(class~., data=training, method="svmRadial", trControl = ctrl,
                     tuneGrid = SVMgrid, preProc = c("scale","YeoJohnson"), verbose=FALSE)

# print(modelSvmRRB)
pred = predict(mod_fit, newdata=testing)
confusionMatrix(data=pred, testing$class)



## Predict mislabelled 
misLabelled <-  data.frame(PC1 = prComp$x[,1][!colsToKeep], 
                           PC2 = prComp$x[,2][!colsToKeep],
                           class =  clin$V3[!colsToKeep])

sum(predict(mod_fit, newdata = misLabelled) != clin$V3[!colsToKeep])/ nrow(misLabelled)

idx <- which(predict(mod_fit, newdata = misLabelled) == clin$V3[!colsToKeep])




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

###
