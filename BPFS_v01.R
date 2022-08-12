library(caret)
library(tidyverse)
library(e1071)
library(factoextra)

## Data Preprocessing
data_prep <- function(data_with_labels){
  data_temp <- data.frame(data_with_labels[,-1])
  
  ## Create train and test sets
  training.samples_1 <- data_temp$label %>%
    createDataPartition(p = 0.8, list = FALSE)
  
  data_train <- data_temp[training.samples_1,]
  data_test <- data_temp[-training.samples_1,]
  
  return(list(data_temp, data_train, data_test))
}

## BPFS with ONLY FIRST PCA

BPFS_1_PCA <- function(data_train, a, b) {
  ## a - parameter for chosing top PCs
  ## b - parameter for chosing top genes from top PCs
  data_temp_1 <- data_train[, -grep('label', colnames(data_train))][,-1] ## Data without labels
  
  ## First PCA on train data
  pca_1 <- prcomp(data_temp_1, scale = FALSE)
  pca_1.var <- pca_1$sdev^2
  pca_1.var.per <- round(pca_1.var/sum(pca_1.var)*100,1)
  
  ## Get top PCAs gene scores
  gene_score_list <- list(abs(pca_1$rotation[,c(1:a)]))
  gene_score_ranked <- lapply(gene_score_list, function(x) sort(x, decreasing = TRUE))
  top_genes_per_pc <- lapply(gene_score_ranked, function(x) names(x[1:floor(length(x)*b)]))
  gene_selected_pca <- unique(unlist(top_genes_per_pc))
  data_temp_2 <- data_train[,match(gene_selected_pca, colnames(data_train))]
  data_temp_2$label <- factor(data_train$label, levels = c(0,1))
  
  return(data_temp_2)
}

## BPFS with TWO PCAs

BPFS <- function(data_train, a, b, t) {
  ## a - parameter for choosing top PCs
  ## b - parameter for choosing top genes from top PCs
  ## t - parameter for choosing top genes from top PCs in second PCA
  data_temp_1 <- data_train[, -grep('label', colnames(data_train))][,-1] ## Data without labels
  
  ## First PCA on train data
  pca_1 <- prcomp(data_temp_1, scale = FALSE)
  pca_1.var <- pca_1$sdev^2
  pca_1.var.per <- round(pca_1.var/sum(pca_1.var)*100,1)
  
  ## Get top PCAs gene scores
  gene_score_list <- list(abs(pca_1$rotation[,1]))
  if (a > 1) {
    for (i in c(2,a)){
      gene_score_list <- append(gene_score_list, list(abs(pca_1$rotation[,i])))
    }
  } else {
    gene_score_list <- list(abs(pca_1$rotation[,1]))
  }
  gene_score_ranked <- lapply(gene_score_list, function(x) sort(x, decreasing = TRUE))
  top_genes_per_pc <- lapply(gene_score_ranked, function(x) names(x[1:floor(length(x)*b)]))
  gene_selected_pca <- unique(unlist(top_genes_per_pc))
  data_temp_2 <- data_train[,match(gene_selected_pca, colnames(data_train))]
  
  ## Transform result dataframe from PCA 1
  data_temp_3 <- as.data.frame(t(data_temp_2))
  
  ## Second PCA on filtered data from PCA 1
  
  pca_2 <- prcomp(data_temp_3, scale = FALSE)
  
  pca_2.var <- pca_2$sdev^2
  pca_2.var.per <- round(pca_2.var/sum(pca_2.var)*100,1)
  
  ## Results from individuals
  
  res.ind <- get_pca_ind(pca_2)
  ind_contrib_list <- list(abs(res.ind$contrib[,1]))
  if (a > 1) {
    for (i in c(2,a)){
      ind_contrib_list <- append(ind_contrib_list, list(abs(res.ind$contrib[,i])))
    }
  } else {
    ind_contrib_list <- list(abs(res.ind$contrib[,1]))
  }
  ind_contrib_ranked <- lapply(ind_contrib_list, function(x) sort(x, decreasing = TRUE)) # contains the contributions (in percentage) of the variables to the principal components
  top_genes_per_pc_2 <- lapply(ind_contrib_ranked, function(x) names(x[1:floor(length(x)*t)]))
  gene_selected_pca_2 <- unique(unlist(top_genes_per_pc_2))
  
  ## Filter with gene_selected_pca_2
  
  data_temp_4 <- data_temp_2[,match(gene_selected_pca_2, colnames(data_temp_2))]
  
  data_temp_4$label <- factor(data_train$label, levels = c(0,1))
  
  return(data_temp_4)
  
}

## SVM
## Factor the target column

SVM_model <- function(data_train_filtered, data_test){
  
  classifier <- svm(formula = label ~ .,
                    data = data_train_filtered,
                    type = "C-classification",
                    kernel="radial"
  )
  
  ## Parameter tuning
  
  # tmodel <- tune(svm,label ~ .,
  #             data = data_train_filtered,
  #             type = "C-classification",
  #             ranges=list(epsilon = seq(0,1,0.1), cost = 2^(2:7)))
  # 
  # mymodel <- tmodel$best.model
  
  pred <- predict(classifier, data_test)
  
  result <- confusionMatrix(table(pred, data_test$label))
  
  return(list(pred, result))
}
