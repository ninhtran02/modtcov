proj <- function(u,v,samples){

  # Computes cov(u,v)/cov(u,u) * u
  # Ninh Tran
  # Last modified
  # 15 May 2023

  return( sum((u*v)/samples)/sum((u*u)/samples) * u )
}

GS_process <- function(mat, samples){

  # Performs Gram-Schmidt process on mat
  # with reference to samples for covariance
  # and variance information
  # Ninh Tran
  # Last modified
  # 15 May 2023

  n.combos <- dim(mat)[2]
  n.rows <- dim(mat)[1]

  for(i in 2:n.rows){
    for(j in (i-1):1){
      mat[i,] <- mat[i,] - proj(mat[j,],mat[i,],samples)
    }
  }
  return(mat)
}


create_contrast.matrix <- function(contrast,samples){

  # Create a contrast matrix where the first row
  # is the main contrast and remaining contrasts
  # form the covariates
  # Ninh Tran
  # Last modified
  # 15 May 2023

  if(length(contrast) <= 2){
    # Not enough columns
    return(rbind(contrast))
  }

  non.zero.contrast.element <- which(contrast != 0)[1]
  mat <- rbind(contrast,diag(length(contrast))[-non.zero.contrast.element,])

  A <- round(GS_process(mat = mat, samples = samples),10)
  A <- A[-1,]
  A <- t(rbind(A))
  colSumsA <- colSums(A)

  zero.col.sum.indices <- which(colSumsA == 0)
  b <- c()
  if(length(zero.col.sum.indices) >= 1){
    for(i in 1:length(zero.col.sum.indices)){
      b.vec <- rep(0,length(colSumsA))
      b.vec[zero.col.sum.indices[i]] <- 1
      b <- rbind(b,b.vec)
    }
  }

  pair.ref <- t(combn(1:length(colSumsA), 2))

  for(i in 1:dim(pair.ref)[1]){
    b.vec <- rep(0,length(colSumsA))
    ref.A <- colSumsA[pair.ref[i,]]
    if(0 %in% ref.A){
      next
    }
    b.vec[pair.ref[i,]] <- c(-ref.A[2]/ref.A[1], 1)
    b <- rbind(b,b.vec)
  }

  m <- round(A %*% t(b),6)
  #m <- m[, qr(m)$pivot[seq_len(qr(m)$rank)]] # Extracts independent columns

  contrast.mat <- rbind(contrast,t(m))
  return(contrast.mat)
}

create_list.of.contrast.matrices <- function(contrast.list, X){

  # Ninh Tran
  # Last modified
  # 15 May 2023

  if(is.numeric(contrast.list)){
    contrast.list <- list(contrast.list)
  }

  list.of.contrast.matrices <- list()
  samples <- colSums(X)
  for(i in 1:length(contrast.list)){
    contrast.mat <- create_contrast.matrix(contrast = contrast.list[[i]], samples = samples)
    list.of.contrast.matrices <- c(list.of.contrast.matrices, list(contrast.mat))
  }

  return(list.of.contrast.matrices)
}

create.siotani.covariates <- function(t_main, covariates_kshirsagar, df){

  # Tranform a matrix of Kshirsagar's t-statistics
  # with Siotani's result
  # Ninh Tran
  # Last modified
  # 15 May 2023

  fact <- ((1 + (t_main^2)/df)^(-1/2))
  covariates_siotani <- covariates_kshirsagar * as.numeric(fact)

  return(covariates_siotani)
}


#' Create moderated t-statistics and covariates
#'
#' @param Y matrix containing log-ratios or log-expression values for a series of arrays, with rows corresponding to genes and columns to samples.
#' @param X design matrix of the microarray experiment, with rows corresponding to samples and columns to coefficients to be estimated. The matrix must have at least 3 columns and the rows must sum to 1.
#' @param cont contrast numeric for the intended hypothesis test.
#' @param lmFit.args list of further arguments for the lmFit function from the limma package.
#' @param eBayes.args list of further arguments for the eBayes function from the limma package.
#'
#' @return A list with the elements
#' \item{mod.t.stat}{numeric of moderated t-statistics.}
#' \item{mod.t.cov}{numeric of moderated t-covariates.}
#' \item{df}{numeric of degrees of freedoms for the null component.}
#' \item{Y} from the from the arguments.
#' \item{X} from the arguments.
#' @export
#'
#' @examples create.t.stat.cov.df(exampleY, exampleX, c(1,-1,0,0,0))

create.t.stat.cov.df <- function(Y,X,cont,
                                 lmFit.args = list(),
                                 eBayes.args = list()){

  # Ninh Tran
  # Last modified
  # 15 May 2023

  if( !all(rowSums(X) == 1) )  stop("rowSums(X) != 1")
  if( ncol(X) < 3 )  stop("ncol(X) < 3")


  lmFit.args <- c(list(object = Y, design = X), lmFit.args)
  fitted.lm <- do.call(what = lmFit, args = lmFit.args)

  cont.matrix <- contrast.matrices <- create_list.of.contrast.matrices(contrast.list = cont,
                                                                X = X)[[1]]

  df <- c()
  stats.mat <- c()
  num.stats <- dim(cont.matrix)[1]
  for(k in 1:num.stats){
    contrast.of.interest <<- cont.matrix[k,]
    contrast.obj <- limma::makeContrasts( contrast.of.interest, levels = colnames(coef(fitted.lm)))
    eBayes.args.current <- c(list(fit = limma::contrasts.fit(fitted.lm, contrast.obj)),eBayes.args)
    ebayesstats <- do.call(what = limma::eBayes, args = eBayes.args.current)
    stats.vec <- ebayesstats$t
    df <- ebayesstats$df.total[1]
    stats.mat <- cbind(stats.mat,stats.vec)
    rm(contrast.obj)
  }

  df <- ebayesstats$df.total

  # Siotani transformation
  covariates_siotani    <- create.siotani.covariates(t_main = stats.mat[,1],
                                            covariates_kshirsagar = stats.mat[,-1],
                                            df = df)


  output.list <- list(mod.t.stat = stats.mat[,1],
                      mod.t.cov = cbind(covariates_siotani),
                      df = df,
                      Y = Y,
                      X = X)
  return(output.list)
}


