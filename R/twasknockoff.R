#' TWAS Knockoff main function
#' @param snpidx A list of numeric vectors, storing the indices of cis-variants for each candidate gene in the risk region.
#' @param ye A list of numeric vectors, storing the gene expression levels for each gene in the eQTL study.
#' @param Xe A list of matrices, storing the cis-genotype matrices for each gene in the eQTL study.
#' @param summarystat Summary statistics for cis-SNPs in the risk region.
#' @param Xp The cis-genotype matrix from GWAS data or reference panel for the risk region.
#' @param removemethod If 'lasso', remove the significant variants (eQTLs) detected in the gene expression prediction model.
#' @param simu If TRUE, print the correlation between predicted gene expression levels and true gene expression levels given as input via yep_true.
#' @param reduced If TRUE, only consider the cis-SNPs within each gene; if FALSE, consider all cis-SNPs in the risk region.
#' @param lambda_r Parameter for the regularization of correlation matrix.
#' @param correlation If 'improved', apply the improved correlation matrix estimation via bootstrap sample.
#' @param nrep Number of bootstrap samples (including the original copy).
#' @param ts Test statistic for knockoffs, selected from {'marginal','susie','lasso','lasso.approx.lambda','squared.zscore'}.
#' @param appr Approximation method for knockoffs.
#' @param yep_true A list of true gene expression levels for GWAS study, only required when simu = TRUE.
#' @param gene_fdr If 'gene': only use genes for FDR control; If 'genesnp': use all genetic elements (including genes and cis-SNPs) for FDR control.
#' @import GhostKnockoff
#' @import glmnet
#' @importFrom stats cov2cor as.dist cor
#' @return A list of three elements, the first element contains test statistics for all genetic elements, the second element contains q-values and knockoff statistics similar to GhostKnockoff, the third element contains SNPs removed from the knockoff procedure (i.e., significant SNPs in gene expression prediction models)
#' @export
TwasKnockoff <- function(snpidx, ye, Xe, summarystat, Xp, removemethod = 'lasso', simu = FALSE, reduced = TRUE, lambda_r = 0.1,
                         correlation = 'improved', nrep = 10, ts = 'lasso',appr = 'sdp', yep_true = NULL, gene_fdr = 'genesnp'){


  standardize <- function(x) {
    (x - mean(x)) / sd(x)
  }
  ye <- lapply(ye, standardize)
  # y <- data$y

  z_scores <- summarystat$b/summarystat$se
  nsample <- mean(summarystat$N)

  Xp <- bedNA(t(Xp))
  Xp <- scale(Xp)
  dim(Xp)


  n <- nrow(Xp)
  p <- ncol(Xp)
  k <- length(ye)


  ye_predicted <- list()

  beta <- matrix(0, nrow=p, ncol=k)

  z_e <- numeric(k)

  errs = list()
  yehats = list()

  removelist <- numeric()


  for(geneid in 1:length(ye)){
    ind_gene <- snpidx[[geneid]]
    X1 <- Xe[[geneid]]
    X1 <- bedNA(t(X1))
    X1 <- scale(X1)
    X2 <- Xp[,ind_gene]

    LD_matrix <- cor(X2)
    #Lasso:
    lasso_cv <- cv.glmnet(X1, ye[[geneid]], alpha = 0.5)
    # Best lambda value
    best_lambda <- lasso_cv$lambda.min
    lasso_mod = glmnet(X1, ye[[geneid]], alpha = 0.5, lambda = best_lambda)

    causal_snp <- ind_gene[which(lasso_mod$beta[, 1] != 0)]
    # print(lasso_mod$beta[which(lasso_mod$beta!=0)])

    if(length(causal_snp)>0){
      removelist <- append(removelist, causal_snp)
      lasso_pred = predict(lasso_mod, newx = X2)

      yehats <- append(yehats, list(predict(lasso_mod, newx = X1)))
      errs <- append(errs, list(ye[[geneid]] - yehats[[geneid]]))

      ye_predicted <- append(ye_predicted, list(lasso_pred))

      if(simu == TRUE & !is.null(yep_true)){
        print(cor(as.numeric(lasso_pred), yep_true[[geneid]]))
      }

      beta[ind_gene,geneid] <- as.numeric(lasso_mod$beta)
    }

    if(length(causal_snp)==0){
      print(paste("gene",geneid,"no significant variant in lasso"))
      ridge_cv <- cv.glmnet(X1, ye[[geneid]], alpha = 0)
      # Best lambda value
      best_lambda <- ridge_cv$lambda.min
      ridge_mod = glmnet(X1, ye[[geneid]], alpha = 0, lambda = best_lambda)
      ridge_pred = predict(ridge_mod, newx = X2)

      yehats <- append(yehats, list(predict(ridge_mod, newx = X1)))
      errs <- append(errs, list(ye[[geneid]] - yehats[[geneid]]))

      if(simu == TRUE & !is.null(yep_true)){
        print(cor(as.numeric(ridge_pred), yep_true[[geneid]]))
      }

      ye_predicted <- append(ye_predicted, list(ridge_pred))
      beta[ind_gene,geneid] <- as.numeric(ridge_mod$beta)
    }


    z_e[geneid] <- t(beta[ind_gene,geneid])%*%z_scores[ind_gene]/sqrt(t(beta[ind_gene,geneid])%*%LD_matrix%*%beta[ind_gene,geneid])
    # z_e[geneid] <- t(beta[[geneid]])%*%z_scores[ind_gene]/sqrt(t(beta[[geneid]])%*%LD_matrix%*%beta[[geneid]])
  }

  # print(z_e[1:length(ye)])

  removelist0 = removelist

  if(reduced == TRUE){
    allsnp <- unique(unlist(snpidx))
    removelist <- match(removelist, allsnp)
    # print(removelist)
    beta = beta[allsnp,]
    LD <- cor(Xp[,allsnp])
  }else{
    LD <- cor(Xp)
  }


  cov1 <- t(beta)%*%LD%*%beta
  cov2 <- LD%*%beta
  cov1_0 <- cov1
  cov2_0 <- cov2

  block1 <- cbind(cov1_0,t(cov2_0))
  block2 <- cbind(cov2_0,LD)
  LD_new_0 <- stats::cov2cor(rbind(block1,block2))

  LD_new <- LD_new_0


  # generate bootstrap samples:


  if(correlation == 'improved'){

    sds = sapply(errs,sd)
    # yphat = ye_predicted
    for (rep in 2:nrep) {
      #if(rep == nrep){
      print(paste("starting iteration",rep))
      #}
      beta_tmp <- matrix(0, nrow=p, ncol=k)
      # ye_predicted_tmp <- matrix(0,nrow=n,ncol=k)
      ye_predicted_tmp <- list()
      for (geneid in 1:length(ye)){

        ind_gene <- snpidx[[geneid]]
        X1 <- Xe[[geneid]]
        X1 <- bedNA(t(X1))
        X1 <- scale(X1)
        X2 <- Xp[,ind_gene]

        # ye_tem <- yehats[[geneid]] + rnorm(length(yehats[[geneid]])) * sds[geneid]
        ye_tem = yehats[[geneid]] + abs(errs[[geneid]]) * rnorm(length(yehats[[geneid]]))

        #Lasso:
        lasso_cv <- cv.glmnet(X1, ye_tem, alpha = 0.5)
        # Best lambda value
        best_lambda <- lasso_cv$lambda.min
        lasso_mod = glmnet(X1, ye_tem, alpha = 0.5, lambda = best_lambda)
        causal_snp <- which(lasso_mod$beta[, 1] != 0)


        if(length(causal_snp)>0){
          lasso_pred = predict(lasso_mod, newx = X2)

          ye_predicted_tmp <- append(ye_predicted_tmp, list(lasso_pred))
          beta_tmp[ind_gene,geneid] <- as.numeric(lasso_mod$beta)

        }

        #### check point

        if(length(causal_snp)==0){
          # print("no significant variant in lasso")
          ridge_cv <- cv.glmnet(X1, ye_tem, alpha = 0)
          # Best lambda value
          best_lambda <- ridge_cv$lambda.min
          ridge_mod = glmnet(X1, ye_tem, alpha = 0, lambda = best_lambda)
          ridge_pred = predict(ridge_mod, newx = X2)


          ye_predicted_tmp <- append(ye_predicted_tmp, list(ridge_pred))
          beta_tmp[ind_gene,geneid] <- as.numeric(ridge_mod$beta)
        }
      }

      if(reduced == TRUE){
        allsnp <- unique(unlist(snpidx))
        beta_tmp = beta_tmp[allsnp,]
      }

      cov1 <- (cov1*(rep-1)+t(beta_tmp)%*%LD%*%beta_tmp)/rep
      cov2 <- (cov2*(rep-1)+LD%*%beta_tmp)/rep
    }

    block1 <- cbind(cov1,t(cov2))
    block2 <- cbind(cov2,LD)
    LD_new <- stats::cov2cor(rbind(block1,block2))

  }

  if(reduced == TRUE){
    allsnp <- unique(unlist(snpidx))
    z_scores = z_scores[allsnp]
  }

  print(z_e)

  z_all <- append(z_e, z_scores)


  if(removemethod=="lasso"){
    index <- unique(removelist)
    index_all <- index+k

    if(length(index_all)>0){
      z_all <- z_all[-index_all]
      LD_new <- LD_new[-index_all,-index_all]
    }

  }

  if(length(z_all) == 1){
    print("only one genetic element left in the model")
    return(NA)
  }

  lambda <- lambda_r
  print(paste0("lambda",lambda))
  LD_modi <- (1-lambda)*LD_new + lambda*Matrix::diag(nrow(LD_new))
  n.study <- nrow(Xp)
  zscores <- matrix(as.vector(z_all), ncol=1)

  LD_modi <- (LD_modi + t(LD_modi)) / 2

  if(ts == 'squared.zscore'){
    fit.prelim <- GhostKnockoff:::GhostKnockoff.prelim(LD_modi, M=1, method=appr, max.size=100)
    GK.stat <- GhostKnockoff:::GhostKnockoff.fit(zscores, n.study, fit.prelim, gamma=1, weight.study=NULL)
    if(gene_fdr == 'gene'){
      GK.filter<-GhostKnockoff.filter(GK.stat$T_0[1:k,], GK.stat$T_k[1:k,])
    }else{
      GK.filter <- GhostKnockoff:::GhostKnockoff.filter(GK.stat$T_0,GK.stat$T_k)
    }
  }else{
    fit.prelim <- GhostKnockoff.prelim(LD_modi, M=1, method=appr, max.size=100)
    print(paste("s:",Matrix::diag(fit.prelim$diag_s)[1:5]))
    GK.stat<-GhostKnockoff.fit(zscores, n.study, fit.prelim, method=ts, type='fdr', M.fwer=1)
    if(gene_fdr == 'gene'){
      GK.filter<-GhostKnockoff.filter(GK.stat$T_0[[1]][1:k],GK.stat$T_k[[1]][1:k,])
    }else{
      GK.filter<-GhostKnockoff.filter(GK.stat$T_0[[1]],GK.stat$T_k[[1]])
    }
  }

  data <- list(GK.stat = GK.stat, GK.filter = GK.filter, removed_list = unique(removelist0))

  return(data)

}

#' Impute missingness
#' @param bed1 genotype matrix
#' @return an imputed matrix
bedNA <- function(bed1){
  for(j in 1:ncol(bed1)){
    temp <- bed1[,j]
    temp[is.na(temp)] <- mean(temp,na.rm=T)
    bed1[,j] <- temp
    #print(j)
  }
  return(bed1)
}






#' TWASKnockoff with multiple knockoffs
#' @param M number of knockoff copies


TwasKnockoff_multi <- function(snpidx, ye, Xe, summarystat, Xp, removemethod = 'lasso', simu = TRUE, reduced = TRUE, lambda_r = 0.1, 
                         correlation = 'improved', nrep = 10, ts = 'lasso',appr = 'sdp', yep_true = NULL, gene_fdr = 'genesnp', M = 5){
  
  print(paste("M = ", M))
  
  library(susieR)
  
  if(ts == 'squared.zscore'){
    library(GhostKnockoff)
    library(glmnet)
    
  }else{
    library(GhostKnockoff)
    library(glmnet)
    library(ghostbasil)
    source("/home/xz527/Rcode/twas_fm/GhostKnockoff.R")
  }
  # library(GhostKnockoff)
  
  standardize <- function(x) {
    (x - mean(x)) / sd(x)
  }
  ye <- lapply(ye, standardize)
  z_scores <- summarystat$b/summarystat$se
  nsample <- mean(summarystat$N)

  
  bedNA <- function(bed1){
    for(j in 1:ncol(bed1)){
      temp <- bed1[,j]
      temp[is.na(temp)] <- mean(temp,na.rm=T)
      bed1[,j] <- temp
      #print(j)
    }
    return(bed1)
  }
  
  Xp <- bedNA(t(Xp))
  Xp <- scale(Xp)
  dim(Xp)
  
  n <- nrow(Xp)
  p <- ncol(Xp)
  k <- length(ye)
  
  ye_predicted <- list()
  
  beta <- matrix(0, nrow=p, ncol=k)
  z_e <- numeric(k)
  errs = list()
  yehats = list()
  removelist <- numeric()
  
  for(geneid in 1:length(ye)){
    ind_gene <- snpidx[[geneid]]
    X1 <- Xe[[geneid]]
    X1 <- bedNA(t(X1))
    X1 <- scale(X1)
    X2 <- Xp[,ind_gene]
    
    LD_matrix <- cor(X2)
    #Lasso:
    lasso_cv <- cv.glmnet(X1, ye[[geneid]], alpha = 0.5)
    # Best lambda value
    best_lambda <- lasso_cv$lambda.min
    lasso_mod = glmnet(X1, ye[[geneid]], alpha = 0.5, lambda = best_lambda)
    
    causal_snp <- ind_gene[which(lasso_mod$beta!=0)]
    # print(lasso_mod$beta[which(lasso_mod$beta!=0)])
    
    if(length(causal_snp)>0){
      removelist <- append(removelist, causal_snp)
      lasso_pred = predict(lasso_mod, newx = X2)
      yehats <- append(yehats, list(predict(lasso_mod, newx = X1)))
      errs <- append(errs, list(ye[[geneid]] - yehats[[geneid]]))
      ye_predicted <- append(ye_predicted, list(lasso_pred))
      
      if(simu == TRUE & !is.null(yep_true)){
        print(cor(as.numeric(lasso_pred), yep_true[[geneid]]))
      }
      beta[ind_gene,geneid] <- as.numeric(lasso_mod$beta)
    }

    if(length(causal_snp)==0){
      print(paste("gene",geneid,"no significant variant in lasso"))
      ridge_cv <- cv.glmnet(X1, ye[[geneid]], alpha = 0)
      # Best lambda value
      best_lambda <- ridge_cv$lambda.min
      ridge_mod = glmnet(X1, ye[[geneid]], alpha = 0, lambda = best_lambda)
      ridge_pred = predict(ridge_mod, newx = X2)
      yehats <- append(yehats, list(predict(ridge_mod, newx = X1)))
      errs <- append(errs, list(ye[[geneid]] - yehats[[geneid]]))
      if(simu == TRUE & !is.null(yep_true)){
        print(cor(as.numeric(ridge_pred), yep_true[[geneid]]))
      }
      
      ye_predicted <- append(ye_predicted, list(ridge_pred))
      beta[ind_gene,geneid] <- as.numeric(ridge_mod$beta)
    }
    z_e[geneid] <- t(beta[ind_gene,geneid])%*%z_scores[ind_gene]/sqrt(t(beta[ind_gene,geneid])%*%LD_matrix%*%beta[ind_gene,geneid])
  }
  
  removelist0 = removelist
  
  if(reduced == TRUE){
    allsnp <- unique(unlist(snpidx))
    removelist <- match(removelist, allsnp)
    # print(removelist)
    beta = beta[allsnp,]
    LD <- cor(Xp[,allsnp])
  }else{
    LD <- cor(Xp)
  }

  cov1 <- t(beta)%*%LD%*%beta
  cov2 <- LD%*%beta
  cov1_0 <- cov1
  cov2_0 <- cov2
  
  block1 <- cbind(cov1_0,t(cov2_0))
  block2 <- cbind(cov2_0,LD)
  LD_new_0 <- cov2cor(rbind(block1,block2))
  LD_new <- LD_new_0
  
  # generate bootstrap samples:
  if(correlation == 'improved' & nrep > 1){
    
    sds = sapply(errs,sd)

    for (rep in 2:nrep) {
      print(paste("starting iteration",rep))

      beta_tmp <- matrix(0, nrow=p, ncol=k)
      ye_predicted_tmp <- list()
      for (geneid in 1:length(ye)){
        
        ind_gene <- snpidx[[geneid]]
        X1 <- Xe[[geneid]]
        X1 <- bedNA(t(X1))
        X1 <- scale(X1)
        X2 <- Xp[,ind_gene]
        
        # ye_tem <- yehats[[geneid]] + rnorm(length(yehats[[geneid]])) * sds[geneid]
        ye_tem = yehats[[geneid]] + abs(errs[[geneid]]) * rnorm(length(yehats[[geneid]]))
        
        #Lasso:
        lasso_cv <- cv.glmnet(X1, ye_tem, alpha = 0.5)
        # Best lambda value
        best_lambda <- lasso_cv$lambda.min
        lasso_mod = glmnet(X1, ye_tem, alpha = 0.5, lambda = best_lambda)
        causal_snp <- which(lasso_mod$beta!=0)
        
        if(length(causal_snp)>0){
          lasso_pred = predict(lasso_mod, newx = X2)
          ye_predicted_tmp <- append(ye_predicted_tmp, list(lasso_pred))
          beta_tmp[ind_gene,geneid] <- as.numeric(lasso_mod$beta)  
        }
        
        if(length(causal_snp)==0){
          # print("no significant variant in lasso")
          ridge_cv <- cv.glmnet(X1, ye_tem, alpha = 0)
          # Best lambda value
          best_lambda <- ridge_cv$lambda.min
          ridge_mod = glmnet(X1, ye_tem, alpha = 0, lambda = best_lambda)
          ridge_pred = predict(ridge_mod, newx = X2)     
          ye_predicted_tmp <- append(ye_predicted_tmp, list(ridge_pred))
          beta_tmp[ind_gene,geneid] <- as.numeric(ridge_mod$beta)
        }
      }
      
      if(reduced == TRUE){
        allsnp <- unique(unlist(snpidx))
        beta_tmp = beta_tmp[allsnp,]
      }  
      cov1 <- (cov1*(rep-1)+t(beta_tmp)%*%LD%*%beta_tmp)/rep
      cov2 <- (cov2*(rep-1)+LD%*%beta_tmp)/rep
    }
    
    block1 <- cbind(cov1,t(cov2))
    block2 <- cbind(cov2,LD)
    LD_new <- cov2cor(rbind(block1,block2))  
  }
  
  if(reduced == TRUE){
    allsnp <- unique(unlist(snpidx))
    z_scores = z_scores[allsnp]
  }
  print(z_e)
  
  z_all <- append(z_e, z_scores)
  
  if(removemethod=="lasso"){
    index <- unique(removelist)
    index_all <- index+k
    if(length(index_all)>0){
      z_all <- z_all[-index_all]
      LD_new <- LD_new[-index_all,-index_all]
    } 
  }
  
  if(length(z_all) == 1){
    print("only one genetic element left in the model")
    return(NA)
  }
  
  print(estimate_s_rss(z_all, LD_new, n = nsample))
  
  lambda <- lambda_r
  print(paste0("lambda",lambda))
  LD_modi <- (1-lambda)*LD_new + lambda*diag(nrow(LD_new))
  n.study <- nrow(Xp)
  zscores <- matrix(as.vector(z_all), ncol=1)
  
  LD_modi <- (LD_modi + t(LD_modi)) / 2
  
  if(ts == 'squared.zscore'){
    fit.prelim <- GhostKnockoff:::GhostKnockoff.prelim(LD_modi, M = M, method=appr, max.size=100)
    GK.stat <- GhostKnockoff:::GhostKnockoff.fit(zscores, n.study, fit.prelim, gamma=1, weight.study=NULL)
    if(gene_fdr == 'gene'){
      GK.filter<-GhostKnockoff.filter(GK.stat$T_0[1:k,], GK.stat$T_k[1:k,])
    }else{
      GK.filter <- GhostKnockoff:::GhostKnockoff.filter(GK.stat$T_0,GK.stat$T_k)
    }
  }else{
    fit.prelim <- GhostKnockoff.prelim(LD_modi, M = M, method=appr, max.size=100)
    print(paste("s:",diag(fit.prelim$diag_s)[1:5]))
    GK.stat<-GhostKnockoff.fit(zscores, n.study, fit.prelim, method=ts, type='fdr', M.fwer=50)
    if(gene_fdr == 'gene'){
      GK.filter<-GhostKnockoff.filter(GK.stat$T_0[[1]][1:k],GK.stat$T_k[[1]][1:k,])
    }else{
      GK.filter<-GhostKnockoff.filter(GK.stat$T_0[[1]],GK.stat$T_k[[1]])
    }
  } 
  
  data <- list(GK.stat, GK.filter, removelist0)
  return(data)
  
}










#' TWASKnockoff using publicly available weights
#' @param p by K beta matrix of weights for eQTLs

TwasKnockoff_weight <- function(snpidx, beta, summarystat, Xp, removemethod = 'lasso', simu = TRUE, reduced = TRUE, lambda_r = 0.1, 
                         correlation = 'improved', nrep = 10, ts = 'lasso',appr = 'sdp', yep_true = NULL, gene_fdr = 'genesnp'){
  
  library(susieR)
  
  if(ts == 'squared.zscore'){
    library(GhostKnockoff)
    library(glmnet)
  }else{
    library(GhostKnockoff)
    library(glmnet)
    library(ghostbasil)
    source("/home/xz527/Rcode/twas_fm/GhostKnockoff.R")
  }
  
  standardize <- function(x) {
    (x - mean(x)) / sd(x)
  }

  z_scores <- summarystat$b/summarystat$se
  nsample <- mean(summarystat$N)
  
  bedNA <- function(bed1){
    for(j in 1:ncol(bed1)){
      temp <- bed1[,j]
      temp[is.na(temp)] <- mean(temp,na.rm=T)
      bed1[,j] <- temp
      #print(j)
    }
    return(bed1)
  }
  
  Xp <- bedNA(t(Xp))
  Xp <- scale(Xp)
  dim(Xp)
  
  n <- nrow(Xp)
  p <- ncol(Xp)
  k <- nrow(risk_gene)
  
  LD_matrix <- cor(Xp)
  z_e <- numeric(k)
  removelist <- numeric()
  
  for(geneid in 1:k){
    removelist <- append(removelist, which(beta[, geneid] != 0))
    z_e[geneid] <- t(beta[,geneid])%*%z_scores/sqrt(t(beta[,geneid])%*%LD_matrix%*%beta[,geneid])
   }
  
  removelist <- unique(removelist)
  removelist0 = removelist
  
  if(reduced == TRUE){
    allsnp <- unique(unlist(snpidx))
    allsnp <- sort(unique(c(allsnp, removelist)))
    removelist <- match(removelist, allsnp)
    beta = beta[allsnp,]
    LD <- cor(Xp[,allsnp])
    z_scores = z_scores[allsnp]
  }else{
    LD <- cor(Xp)
  }
  
  
  cov1 <- t(beta)%*%LD%*%beta
  cov2 <- LD%*%beta
  cov1_0 <- cov1
  cov2_0 <- cov2
  
  block1 <- cbind(cov1_0,t(cov2_0))
  block2 <- cbind(cov2_0,LD)
  LD_new_0 <- cov2cor(rbind(block1,block2))
  
  LD_new <- LD_new_0
  
  
  
  print(z_e)
  
  z_all <- append(z_e, z_scores)
  
  
  if(removemethod=="lasso"){
    index <- unique(removelist)
    index_all <- index+k
    
    if(length(index_all)>0){
      z_all <- z_all[-index_all]
      LD_new <- LD_new[-index_all,-index_all]
    }
    
  }
  
  if(length(z_all) == 1){
    print("only one genetic element left in the model")
    return(NA)
  }
  
  print(estimate_s_rss(z_all, LD_new, n = nsample))
  
  lambda <- lambda_r
  print(paste0("lambda",lambda))
  LD_modi <- (1-lambda)*LD_new + lambda*diag(nrow(LD_new))
  n.study <- nrow(Xp)
  zscores <- matrix(as.vector(z_all), ncol=1)
  
  LD_modi <- (LD_modi + t(LD_modi)) / 2
  
  
  if(ts == 'squared.zscore'){
    fit.prelim <- GhostKnockoff:::GhostKnockoff.prelim(LD_modi, M=1, method=appr, max.size=100)
    GK.stat <- GhostKnockoff:::GhostKnockoff.fit(zscores, n.study, fit.prelim, gamma=1, weight.study=NULL)
    if(gene_fdr == 'gene'){
      GK.filter<-GhostKnockoff.filter(GK.stat$T_0[1:k,], GK.stat$T_k[1:k,])
    }else{
      GK.filter <- GhostKnockoff:::GhostKnockoff.filter(GK.stat$T_0,GK.stat$T_k)
    }
  }else{
    fit.prelim <- GhostKnockoff.prelim(LD_modi, M=1, method=appr, max.size=100)
    print(paste("s:",diag(fit.prelim$diag_s)[1:5]))
    GK.stat<-GhostKnockoff.fit(zscores, n.study, fit.prelim, method=ts, type='fdr', M.fwer=1)
    if(gene_fdr == 'gene'){
      GK.filter<-GhostKnockoff.filter(GK.stat$T_0[[1]][1:k],GK.stat$T_k[[1]][1:k,])
    }else{
      GK.filter<-GhostKnockoff.filter(GK.stat$T_0[[1]],GK.stat$T_k[[1]])
    }
  } 
  
  
  data <- list(GK.stat, GK.filter, removelist0)
  return(data)
  
}



