#' @useDynLib perms C_get_log_ML
#' @useDynLib perms C_get_log_ML_bioassay

#' @title get_log_ML
#' @description Computes the log marginal likelihood of the data from the log permanents.
#' Given the computed log permanents log_perms, this function computes the log marginal 
#' likelihood using the formula (2.3) in [1]. It is assumed that there are no repeated trials.
#' If the data contain repeated trials, then the appropriate log binomial factor must be added 
#' to the output of this function.
#' @param log_perms A vector length n containing the computed log permanents, where a zero permanent
#'  is indicated by a NA value.
#' @param n Sample size. 
#' @param debug If \code{TRUE}, debug information is printed. 
#' @return The estimated log marginal likelihood. A NA value is returned if there are no non-zero
#' numbers. 
#' @examples
#' library(perms)
#' set.seed(1996)
#' n = 100
#' t = seq(0, 1, length.out=n)
#' y = c(rep(0, n/2), rep(1, n/2))
#' S = 200
#' X = matrix(runif(n*S),nrow = S, ncol = n)
#'
#' log_perms = get_log_perms(X, t, y, debug = FALSE, parallel = TRUE, num_cores = NULL)
#'
#' num_nonzero_perms = sum(!is.na(log_perms))
#' num_nonzero_perms
#'
#' log_ML = get_log_ML(log_perms, n, FALSE)
#' log_ML
#' @references
#' [1] Christensen, D (2023). Inference for Bayesian nonparametric models with binary response data via permutation counting. Bayesian Analysis, Advance online publication, DOI: 10.1214/22-BA1353.
#' @export
get_log_ML= function(log_perms, n, debug=FALSE){
  
  if(!is.vector(log_perms)){
    stop("log_perms must be a vector")
  }
  S = length(log_perms)
  
  res = .Call(C_get_log_ML, as.numeric(log_perms),  as.integer(n), as.integer(S), as.integer(debug))
  
  return(res)
}


#' @title get_log_ML_bioassay
#' @description Computes the log marginal likelihood of bioassay data from the log permanents.
#' Given the computed log permanents log_perms, this function computes the log marginal likelihood
#' using the formula (2.3) in [1]. It takes care of repeated trials by adding the appropriate 
#' log binomial factor.
#' @param log_perms A vector length n containing the computed log permanents, where a zero permanent
#'  is indicated by a NA value.
#' @param successes A vector of length n contatining the number of successful trials at
#' each level.
#' @param trials A vector of length n containing the number of trials at each level.
#' @param n Sample size. 
#' @param debug If \code{TRUE}, debug information is printed. 
#' @return The estimated log marginal likelihood. A NA value is returned if there are no non-zero
#' numbers. 
#' @examples
#' ## Dirichlet toy model
#' library(perms)
#' set.seed(1996)
#' n = 500
#' num_trials = 10
#' levels = seq(-1, 1, length.out = num_trials)
#' 
#' trials = rep(n %/% num_trials, num_trials)
#' successes = c(10, 26, 10, 20, 20, 19, 29, 24, 31, 33)
#' 
#' S = 200
#' alpha = 1.0
#' 
#' get_X = function(S,n,alpha,seed){
#'   set.seed(seed)
#'   X = matrix(0, nrow = S, ncol = n)
#'   for (s in 1:S) {
#'     X[s,1] = rnorm(1)
#'     for (i in 2:n) {
#'       u = runif(1)
#'       if(u < (alpha/(alpha+i-1))){
#'         X[s,i] = rnorm(1)
#'       }else{
#'         if(i==2){
#'           X[s,i] = X[s,1]
#'         }else{
#'           X[s,i] = sample(X[s, 1:(i-1)],size=1)
#'         }
#'       }
#'       
#'     }
#'     
#'   }
#'   return(X)
#' }
#' 
#' seed = 1996
#' X = get_X(S, n, alpha, seed)
#' log_perms = get_log_perms_bioassay(X, levels, successes, trials,
#'            debug=FALSE,parallel = FALSE)
#' log_ml = get_log_ML_bioassay(log_perms, successes, trials, n)
#' 
#' proportion = sum(!is.na(log_perms)) / S*100
#' 
#' proportion 
#' log_ml
#' @references
#' [1] Christensen, D (2023). Inference for Bayesian nonparametric models with binary response data via permutation counting. Bayesian Analysis, Advance online publication, DOI: 10.1214/22-BA1353.
#' @export
get_log_ML_bioassay= function(log_perms, successes, trials, n, debug=FALSE){
  
  if(!is.vector(log_perms)){
    stop("log_perms must be a vector")
  }
  if( !is.vector(successes) || !is.vector(trials)){
    stop("successes and trials must all be vectors")
  }
  
  S = length(log_perms)
  num_trials = length(successes)
  
  if(length(trials)!= num_trials){
    stop("trials must have same length as successes")
  }
  
  if(sum(trials)!= n){
    stop("The sum of trials must be equal to n")
  }
  res = .Call(C_get_log_ML_bioassay, as.numeric(log_perms),  as.integer(successes), as.integer(trials), 
              as.integer(n), as.integer(num_trials),as.integer(S), as.integer(debug))
  
  return(res)
}





