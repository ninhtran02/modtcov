#' Create mimicked local false discovery rates
#'
#' @param mod.t.stat numeric of moderated t-statistics.
#' @param prob_mat matrix of the conditional probabilities.
#' @param param numeric of fitted mu and 1/tau parameters.
#' @param df numeric of degrees of freedoms for the null component.
#' @param alt.df numeric of degrees of freedoms for the non-null components.
#'
#' @return numeric of mimicked local false discovery rates.
#' @export
#'
#' @examples mod.t.cov.df <- create.t.stat.cov.df(exampleY, exampleX, c(1,-1,0,0,0))
#' @examples fit <- fit.t.model(mod.t.cov.df)
#' @examples mlfdr <- mlfdr_fn_tmm(mod.t.stat = mod.t.cov.df$mod.t.stat, prob_mat = fit$prob_mat, param = fit$param, df = mod.t.cov.df$df_vec, alt.df = mod.t.cov.df$df_vec)
mlfdr_fn_tmm <- function(mod.t.stat, prob_mat,
                         param, df, alt.df){

  # Fits mimicked lfdr statistics.
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

  # Compute the pdfs
  mlfdr_main_stat_vec  <- pr_mat * (1/sigma_mat) * dt(x = (t_mat - mu_mat)/sigma_mat, df = df_mat)
  # Compute the lfdrs
  mlfdr_main_stat <-  mlfdr_main_stat_vec[,n.comp+1]/rowSums(mlfdr_main_stat_vec)

  return(mlfdr_main_stat)
}
