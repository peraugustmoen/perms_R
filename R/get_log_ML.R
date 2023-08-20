#' @useDynLib perms C_get_log_ML
#' @useDynLib perms C_get_log_ML_bioassay

#' @title get_log_ML
#' @description Computes the log marginal likelihood of the data from the log permanents.
#' Given the computed log permanents logperms, this function computes the log marginal 
#' likelihood using the formula (2.3) in [1]. It is assumed that there are no repeated trials.
#' If the data contain repeated trials, then the appropriate log binomial factor must be added 
#' to the output of this function.
#' @param logperms A vector length n containing the computed log permanents, where a zero permanent
#'  is indicated by a -1.
#' @param n Sample size. 
#' @param S Number of samples from the data model. That is, the number of iterations in the estimator.
#' @param debug If \code{TRUE}, debug information is printed. 
#' @return The estimated log marginal likelihood.
#' @examples
#' @references
#' [1] Christensen, D (2023). Inference for Bayesian nonparametric models with binary response data via permutation counting. Bayesian Analysis, Advance online publication, DOI: 10.1214/22-BA1353.
#' @export
get_log_ML= function(logperms, n, S, debug=FALSE){
  
  if(!is.vector(logperms)){
    stop("logperms must be a vector")
  }
  if(length(logperms)!= S){
    stop("logperms must have length S")
  }
  
  res = .Call(C_get_log_ML, as.numeric(logperms),  as.integer(n), as.integer(S), as.integer(debug))
  
  return(res)
}


#' @title get_log_ML_bioassay
#' @description Computes the log marginal likelihood of bioassay data from the log permanents.
#' Given the computed log permanents logperms, this function computes the log marginal likelihood
#' using the formula (2.3) in [1]. It takes care of repeated trials by adding the appropriate 
#' log binomial factor.
#' @param logperms A vector length n containing the computed log permanents, where a zero permanent
#'  is indicated by a -1.
#' @param successes A vector of length n contatining the number of successful trials at
#' each level.
#' @param trials A vector of length n containing the number of trials at each level.
#' @param n Sample size. 
#' @param num_trials Number of different trials.  
#' @param S Number of samples from the data model. That is, the number of iterations in the estimator.
#' @param debug If \code{TRUE}, debug information is printed. 
#' @return The estimated log marginal likelihood.
#' @examples
#' @references
#' [1] Christensen, D (2023). Inference for Bayesian nonparametric models with binary response data via permutation counting. Bayesian Analysis, Advance online publication, DOI: 10.1214/22-BA1353.
#' @export
get_log_ML_bioassay= function(logperms, successes, trials, n, num_trials, S, debug=FALSE){
  
  if(!is.vector(logperms)){
    stop("logperms must be a vector")
  }
  if(length(logperms)!= S){
    stop("logperms must have length S")
  }
  if(length(successes)!= num_trials){
    stop("successes must have length num_trials")
  }
  if(length(trials)!= num_trials){
    stop("trials must have length num_trials")
  }
  if(sum(trials)!= n){
    stop("the sum of trials must be equal to n")
  }
  res = .Call(C_get_log_ML_bioassay, as.numeric(logperms),  as.integer(successes), as.integer(trials), 
              as.integer(n), as.integer(num_trials),as.integer(S), as.integer(debug))
  
  return(res)
}





