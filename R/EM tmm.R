## Hmat function
form_Hmat_tmm <- function(param,
                          prob_mat,
                          mod.t.stat,
                          df,
                          alt.df){

  # The "E" step of the EM algorithm.
  # Ninh Tran
  # Last modified
  # 15 May 2023

  # Number of alternative density components
  n.comp <- length(param)/2

  # Create a matrix of mu values
  mu_mat    <- cbind(matrix(param[1:n.comp],
                            nrow = length(mod.t.stat),
                            ncol = n.comp, byrow = TRUE),0)
  # Create a matrix of sigma values
  sigma_mat <-cbind(matrix(param[(n.comp + 1):(2*n.comp)],
                           nrow = length(mod.t.stat),
                           ncol = n.comp, byrow = TRUE),1)
  # Create a matrix of weights
  pr_mat <- cbind(prob_mat,1-rowSums(prob_mat))
  # Create a matrix of moderate t-statistics
  t_mat <- matrix(mod.t.stat, nrow = length(mod.t.stat), ncol = n.comp+1)
  # Create a matrix of degrees of freedoms
  df_mat <- cbind(matrix(alt.df, ncol = n.comp, nrow = length(mod.t.stat)), df)
  # Compute the computational probabilities
  density_mat <- pr_mat * (1/sigma_mat) * dt(x = (t_mat - mu_mat)/sigma_mat, df = df_mat)
  Hmat <- density_mat/rowSums(density_mat)

  return(Hmat)
}

## Objective functions
neg_Q_fn_tmm <- function(param,
                         Hmat,
                         mod.t.stat,
                         alt.df){

  # The objective function for the "M" step of the EM algorithm
  # Ninh Tran
  # Last modified
  # 15 May 2023

  # Number of alternative density components
  n.comp <- length(param)/2

  # Create a matrix of mu values
  mu_mat    <- matrix(param[1:n.comp],
                      nrow = length(mod.t.stat),
                      ncol = n.comp, byrow = TRUE)
  # Create a matrix of sigma values
  sigma_mat <- matrix(param[(n.comp + 1):(2*n.comp)],
                      nrow = length(mod.t.stat),
                      ncol = n.comp, byrow = TRUE)
  # Create a matrix of moderate t-statistics
  t_mat <- matrix(mod.t.stat, nrow = length(mod.t.stat), ncol = n.comp)
  # Create a matrix of degrees of freedoms
  df_mat <- matrix(alt.df, ncol = n.comp, nrow = length(mod.t.stat))
  # Compute the negative incomplete loglikelihood
  loglike_mat <- -log(sigma_mat) + dt(x = (t_mat - mu_mat)/sigma_mat, df = df_mat, log = TRUE)
  loglike_mat <- loglike_mat*Hmat[,-(n.comp + 1)]
  Q <- -sum(loglike_mat)

  return(Q)
}

# Gradient functions
neg_Q_gn_tmm <- function(param,
                         Hmat,
                         mod.t.stat,
                         alt.df){

  # The derivative of the objective function in the "M" step.
  # Ninh Tran
  # Last modified
  # 15 May 2023

  # Number of alternative density components
  n.comp <- length(param)/2

  # Create a matrix of mu values
  mu_mat    <- matrix(param[1:n.comp],
                      nrow = length(mod.t.stat),
                      ncol = n.comp, byrow = TRUE)
  # Create a matrix of sigma values
  sigma_mat <- matrix(param[(n.comp + 1):(2*n.comp)],
                      nrow = length(mod.t.stat),
                      ncol = n.comp, byrow = TRUE)
  # Create a matrix of moderate t-statistics
  t_mat <- matrix(mod.t.stat, nrow = length(mod.t.stat), ncol = n.comp)
  # Create a matrix of degrees of freedoms
  df_mat <- matrix(alt.df, ncol = n.comp, nrow = length(mod.t.stat))

  # Compute derivatives
  subkernel <- (t_mat - mu_mat)/sigma_mat
  kernel <- (1 + (subkernel^2)/df_mat )
  beta_kernel <- ((df_mat+1)*(subkernel)/(sigma_mat * df_mat * kernel))
  sigma_kernel <- ((df_mat+1)*(subkernel^2)/(sigma_mat * df_mat * kernel)) - 1/sigma_mat
  beta_grad  <- Hmat[,-(n.comp+1)] * (beta_kernel)
  sigma_grad <- Hmat[,-(n.comp+1)] * (sigma_kernel)
  Q_grad <- -colSums(cbind(beta_grad,sigma_grad))

  return(Q_grad)
}



# Parameter fitting functions
EM_tmm = function(mod.t.stat,
                  param, prob_mat,
                  df, alt.df, tol = 10e-6, mod.t.cov,
                  n.iter = NULL, max.iter = 25, n.nodes = NULL){

  # Executes the EM algorithm.
  # Ninh Tran
  # Last modified
  # 16 May 2023

  if(is.null(n.iter)) n.iter <- min(4*(ncol(mod.t.cov)+2), max.iter)

  # First training
  output <- list(param = param, prob_mat = prob_mat)
  output <- EM_fix_point_fn_tmm(param = output$param, prob_mat = output$prob_mat,
                                mod.t.stat = mod.t.stat, mod.t.cov =  mod.t.cov,
                                df = df, alt.df = alt.df, tol = tol, n.nodes = n.nodes)
  #ll_prev <- likelihood_tmm(param = output$param,
  #                          prob_mat = output$prob_mat,
  #                          mod.t.stat = mod.t.stat,
  #                          df = df,
  #                          alt.df = alt.df)
  # Further training
  for(i in 1:(n.iter-1) ){
    output <- EM_fix_point_fn_tmm(param = output$param, prob_mat = output$prob_mat,
                                  mod.t.stat = mod.t.stat, mod.t.cov =  mod.t.cov,
                                  df = df, alt.df = alt.df, tol = tol,
                                  n.nodes = n.nodes)
    #ll_cur <- likelihood_tmm(param = output$param,
    #                         prob_mat = output$prob_mat,
    #                         mod.t.stat = mod.t.stat,
    #                         df = df,
    #                         alt.df = alt.df)

    #if( ll_cur - ll_prev > 0 & abs(ll_cur - ll_prev) < tol ){
    #  break
    #}
    #ll_prev <- ll_cur
  }

  return(output)
}

EM_fix_point_fn_tmm = function(param, prob_mat, mod.t.stat, mod.t.cov,
                               df, alt.df, tol, n.nodes = NULL){

  # Executes the "E" and "M" steps of the EM algorithm.
  # Ninh Tran
  # Last modified
  # 15 May 2023

  k <- ncol(mod.t.cov)
  if(is.null(n.nodes)) {n.nodes <- max(2*ncol(mod.t.cov),4); n.nodes <- min(n.nodes,25)}
  n.comp <- length(param)/2

  Hmat <- form_Hmat_tmm(param = param, prob_mat = prob_mat,
                        mod.t.stat = mod.t.stat,
                        df = df, alt.df)

  mod.t.cov.scaled <- scale(x = cbind(mod.t.cov), center = FALSE, scale = FALSE)
  colnames(Hmat) <- paste("prob",1:(n.comp+1),sep = "")
  colnames(mod.t.cov.scaled) <- paste("x",1:k,sep = "")
  data <- as.data.frame(cbind(Hmat,mod.t.cov.scaled))
  formula <- paste(paste(colnames(Hmat),collapse = "+"), paste(colnames(mod.t.cov.scaled),collapse = "+"),sep = "~")

  nnobj <- nnet::nnet(x = mod.t.cov, y = Hmat, size = n.nodes,
                      softmax = TRUE, maxit = 200, trace = FALSE)
  pr_mat <- predict(object = nnobj, newdata = mod.t.cov)

  param <- optim(par = param[1:(2*n.comp)],
                 fn = neg_Q_fn_tmm,
                 gr = neg_Q_gn_tmm,
                 mod.t.stat = mod.t.stat , Hmat = Hmat,
                 alt.df = alt.df,
                 method = "L-BFGS-B", lower = c(rep(-Inf,n.comp),rep(10e-4,n.comp)),
                 control = list(factr = tol))$par

  output <- list(param = param,
                 prob_mat = pr_mat[,1:n.comp],
                 nnobj = nnobj,
                 mod.t.stat = mod.t.stat,
                 mod.t.cov = mod.t.cov,
                 df = df,
                 alt.df = alt.df)
  return(output)
}


likelihood_tmm <- function(param,
                           prob_mat,
                           mod.t.stat,
                           df,
                           alt.df){

  # Computes the loglikelihood of the tmm.
  # Ninh Tran
  # Last modified
  # 15 May 2023

  # Number of alternative density components
  n.comp <- length(param)/2

  # Create a matrix of mu values
  mu_mat    <- cbind(matrix(param[1:n.comp],
                            nrow = length(mod.t.stat),
                            ncol = n.comp, byrow = TRUE),0)
  # Create a matrix of sigma values
  sigma_mat <-cbind(matrix(param[(n.comp + 1):(2*n.comp)],
                           nrow = length(mod.t.stat),
                           ncol = n.comp, byrow = TRUE),1)
  # Create a matrix of weights
  pr_mat <- cbind(prob_mat,1-rowSums(prob_mat))
  # Create a matrix of moderate t-statistics
  t_mat <- matrix(mod.t.stat, nrow = length(mod.t.stat), ncol = n.comp+1)
  # Create a matrix of degrees of freedoms
  df_mat <- cbind(matrix(alt.df, ncol = n.comp, nrow = length(mod.t.stat)), df)
  # Compute the computational probabilities
  density_mat <- pr_mat * (1/sigma_mat) * dt(x = (t_mat - mu_mat)/sigma_mat, df = df_mat)
  loglik <- sum(log(rowSums(density_mat)))

  return(loglik)
}

#' Fit a neural network t-mixture model
#'
#' @param mod.t.cov.df list returned from \link[modtcov]{create.t.stat.cov.df}.
#' @param tol tolerance control for the "L-BFGS-B" method.
#' @param alt.df numeric of the degrees of freedom for the non-null t-mixture model components. Taken to be the same as the null component by default.
#' @param n.iter number of iterations for the EM algorithm.
#' @param max.iter maximum number of iterations for the EM algorithm.
#' @param n.nodes number of nodes in the hidden layer of the neural network.
#'
#' @return A list with the elements
#' \item{param}{numeric of fitted mu parameters and 1/tau parameters.}
#' \item{prob_mat}{numeric of fitted conditional probabilities.}
#' \item{nnobj}{object from \link[nnet]{nnet}.}
#' \item{mod.t.stat}{numeric of moderated t-statistics.}
#' \item{mod.t.cov}{numeric of moderated t-covariates.}
#' \item{df}{numeric of degrees of freedoms for the null component.}
#' \item{alt.df}{numeric of degrees of freedoms for the non-null components.}
#' @export
#'
#' @examples mod.t.cov.df <- create.t.stat.cov.df(exampleY, exampleX, c(1,-1,0,0,0))
#' @examples fit <- fit.t.model(mod.t.cov.df)
fit.t.model <- function(mod.t.cov.df,
                        tol = 1e-7,
                        alt.df = NULL,
                        n.iter = NULL,
                        max.iter = 25,
                        n.nodes = NULL){

  # Fits tmm with unsupervised initial parameter selections.
  # Ninh Tran
  # Last modified
  # 20 May 2023


  # Get moderate t statistics and covariates
  mod.t.cov <- mod.t.cov.df$mod.t.cov
  mod.t.stat <- mod.t.cov.df$mod.t.stat
  df <- mod.t.cov.df$df
  # Get the d.o.f. for the alternative components
  if(is.null(alt.df)) alt.df <- mod.t.cov.df$df

  # Number of genes
  m <- nrow(mod.t.cov)
  # Number of treatments
  p <- ncol(X)
  # Temporary number of t-mixture components
  Rtemp <- p*3
  # Guess the null proportion
  p.value <- 2-2*pt(q = abs(mod.t.stat), df = df)
  pi0_guess <- limma::propTrueNull(p.value, method="lfdr", nbins=20)*0.6

  # Compute the number of t-mixture components
  clusterobj <- kmeans(x = mod.t.stat, centers = Rtemp)
  cluster.center.order <- order(clusterobj$centers^2)
  null.center.index <- which(cumsum(clusterobj$size[cluster.center.order])/m > pi0_guess)[1]
  if( Rtemp - (null.center.index-1) < 2  ) null.center.index <- Rtemp-1
  non.null.centers <- clusterobj$centers[cluster.center.order[-(1:(null.center.index-1))]]
  R <- length(non.null.centers) + 1

  # Compute the initial parameters
  mu_vec <- non.null.centers
  sigma_vec <- rep(sd(mod.t.stat),length(mu_vec))
  init_prob <- clusterobj$size[cluster.center.order[-(1:(null.center.index-1))]]/m
  prob_mat <- matrix(data = init_prob, nrow = m, ncol = length(init_prob), byrow = TRUE)

  param <- c(mu_vec, sigma_vec)
  n.comp <- R-1

  fit <- EM_tmm(mod.t.stat = mod.t.stat,
                      param = param, prob_mat = prob_mat, df = df,
                      alt.df = alt.df, tol = tol, mod.t.cov = mod.t.cov,
                      n.iter = n.iter, max.iter = max.iter)
  return(fit)
}

#' Generate moderated t-statistics from the neural network t-mixture model
#'
#' @param fit list returned by \link[modtcov]{fit.t.model}.
#'
#' @return
#' @export
#'
#' @examples mod.t.cov.df <- create.t.stat.cov.df(exampleY, exampleX, c(1,-1,0,0,0))
#' @examples fit <- fit.t.model(mod.t.cov.df)
#' @examples generate_modt(fit)
generate_modt <- function(fit){

  # Generate moderated t-statistics using the fitted tmm.
  # Ninh Tran
  # Last modified
  # 15 May 2023

  gt_vec <- c()
  prob_mat <- predict(object = fit$nnobj, newdata = as.data.frame(fit$mod.t.cov))
  for(i in 1:dim(fit$mod.t.cov)[1]){
    # Choose a component
    r <- sample(x = 1:dim(prob_mat)[2], size = 1, replace = TRUE, prob = prob_mat[i,])
    # Select the parameters of that component
    mu <- c(fit$param[1:(length(fit$param)/2)],0)[r]
    sigma <- c(fit$param[(length(fit$param)/2 + 1):length(fit$param)],1)[r]
    df.t <- c(rep(fit$alt.df[i], (length(fit$param)/2)), fit$df[i] )[r]
    # Generate a moderated t-statistic
    gt <- sigma*rt(n = 1, df = df.t) + mu
    gt_vec <- c(gt_vec, gt)
  }
  return(gt_vec)
}
