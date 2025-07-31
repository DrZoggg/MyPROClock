





fitness.AUC_zg <- function (combine.dysreg, dysreg, dd5, dd6)
{
  library(dplyr);library(glmnet)
  regs <- dysreg[combine.dysreg != 0, ]
  genes <- union(regs[, 1], regs[, 2])

  dd7 <- list()
  for(j in names(dd5)){
    ExpData_c <- dd5[[j]]
    ClinData_c <- dd6[[j]]
    pheno.data_c <- ClinData_c[,c("id", "State")]
    pheno.data_c <- pheno.data_c[!is.na(pheno.data_c$State),]
    colnames(pheno.data_c) <- c('sample',"clin.factor")
    exp <- ExpData_c[rownames(ExpData_c) %in% genes, ]
    exp <- as.data.frame(t(exp), stringsAsFactors = F)
    exp$sample <- rownames(exp)
    colnames(exp) <- sub("-", "_", colnames(exp))
    pheno.data_c <- pheno.data_c[ !is.na(pheno.data_c$clin.factor),  ]
    data <- merge(exp, pheno.data_c, by = "sample")
    data <- data[, -1]
    dd7[[j]] <- data
}


  cv.auc <- data.frame()
  for(i in names(dd7)){
    data <- dd7$GSE59867
    train.data <- data
    test.data <- dd7[[i]]
    genes <- setdiff(colnames(test.data), c("clin.factor"))
    x <- train.data[, colnames(train.data) %in% genes] %>%as.matrix()
    y <- train.data$clin.factor
    y <- ifelse(y%in%'Deterioration',1,0)
    if(ncol(x)>=2){

      fit <- glmnet(x, y, family = "gaussian",nlambda = 100, alpha = 1)
 
      cvfit <- cv.glmnet(x, y,
                         family = 'gaussian',
                         alpha = 1,
                         nfolds = 3,
                         type.measure="mse")
      coef.min = coef(cvfit, s = "lambda.min") 
      active.min = which(as.numeric(coef.min)!=0)
      geneids <- colnames(x)[active.min[-1]-1]

      score <- predict(fit,newx = as.matrix(test.data[,colnames(test.data)%in%
                                                        genes]),
                       s = cvfit$lambda.min,
                       type = 'response')%>%as.data.frame()
      score$clin.factor <- test.data$clin.factor
      score$clin.factor <- ifelse(score$clin.factor%in%'Deterioration',1,0)
      auc.ci <- pROC::roc(score$clin.factor,score$s1)$auc%>%as.numeric()%>%suppressMessages();auc.ci
    }else{
      score <- data.frame(s1 = test.data[,1])
      score$clin.factor <- test.data$clin.factor
      score$clin.factor <- ifelse(score$clin.factor%in%'Deterioration',1,0)
      auc.ci <- pROC::roc(score$clin.factor,score$s1)$auc%>%as.numeric()%>%suppressMessages();auc.ci
    }

    res.i <- data.frame(batch = i, fitness = auc.ci, stringsAsFactors = F)
    cv.auc <- rbind(cv.auc,res.i)
}

cv.auc$times <- c(1:nrow(cv.auc))
cv.auc <- mean(cv.auc$fitness)

return(cv.auc)

}


