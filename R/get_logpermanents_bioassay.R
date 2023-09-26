#' @useDynLib perms C_get_log_perms_bioassay
#' @import doParallel
#' @import parallel
#' @import foreach
#' @title get_log_perms_bioassay
#' @description Computes log permanents associated with simulated latent variables X with bioassay data.
#' Each row of the matrix X contains a random sample of size n from the data model. 
#' The observed data are represented as (levels, successes, trials), where levels are the different
#'  levels at which trials were conducted, successes is the vector of the number of successes per level, 
#'  and trials is the vector of the total number of trials per level. The function returns a vector of
#'  log permanents corresponding to each sample. Note that n must be equal to the sum of the entries of trials.
#' @param X A matrix of dimension S x n, in which each row contains a sample from
#' the data model.
#' @param levels A vector containing the levels at which trials were conducted.
#' @param successes A vector containing the number of successful trials at each level. 
#' @param trials A vector containing the number of trials at each level.
#' @param debug If \code{TRUE}, debug information is printed. 
#' @param parallel If \code{TRUE}, computation is run on several cores
#' @param num_cores (Optional) Specifies the number of cores to use if \code{parallel = TRUE}
#' @return Numpy array of log permanents, each element associated to the corresponding row in X.
#' A zero valued permanent is indicated by a NA value.
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
#' proportion = sum(!is.na(log_perms)) / S*100
#' 
#' proportion 
#' @references
#' [1] Christensen, D (2023). Inference for Bayesian nonparametric models with binary response data via permutation counting. Bayesian Analysis, Advance online publication, DOI: 10.1214/22-BA1353.
#' @export
get_log_perms_bioassay= function(X, levels, successes, trials, debug=FALSE, parallel = TRUE, num_cores = NULL){
  
  # X is n times T
  n = 1
  num_trials = 1
  S = 1
  if(is.vector(X)){
    
    if(debug){
      cat("X is input as vector\n")
    }
    n = length(X)
    S = 1
    
  }else if(is.matrix(X)){
    if(debug){
      cat("X is input as matrix\n")
    }
    S = dim(X)[1]
    n = dim(X)[2]
  }else{
    stop("Input X must be an S x n matrix, or a vector of length n if S = 1")
  }
  
  if(!is.vector(levels) || !is.vector(successes) || !is.vector(trials)){
    stop("levels, successes, trials must all be vectors")
  }
  
  num_trials = length(levels)
  
  if(length(trials)!= num_trials){
    stop("trials must have the same length as levels")
  }
  if(length(successes)!= num_trials){
    stop("successes must have the same length as levels")
  }
  if(sum(trials)!= n){
    stop("The sum of all elements of trials must be equal to n")
  }
  
  if(is.vector(X)){
    res = .Call(C_get_log_perms_bioassay, as.numeric(X[]), as.numeric(levels), as.integer(successes), as.integer(trials), as.integer(n), as.integer(num_trials),as.integer(S), as.integer(debug))
  }
  else{
    if(!parallel){
      res = .Call(C_get_log_perms_bioassay, as.numeric(t(X[])),as.numeric(levels), as.integer(successes), as.integer(trials), as.integer(n), as.integer(num_trials),as.integer(S), as.integer(debug))
    }else{

      if(is.null(num_cores)){
        num_cores <- detectCores()
        if(debug){
          cat(sprintf("Registered %d cores\n", as.integer(num_cores)))
        }
      }
      
      Sdiv = floor(S / num_cores)
      rest = S - Sdiv * num_cores
      starts = 1:num_cores
      stops = 1:num_cores
      
      if(rest==0){
        starts = (1:num_cores - 1)*Sdiv +1
        stops = starts + Sdiv-1
      }
      else{
        starts[1:rest] = (1:rest - 1)*(Sdiv+1) +1
        starts[(rest+1):num_cores] = ((rest+1):num_cores - 1)*(Sdiv) +rest +1
        stops[1:rest] = starts[1:rest]+Sdiv
        stops[(rest+1):num_cores] = starts[(rest+1):num_cores]+Sdiv-1
      }
      
      registerDoParallel(cores=num_cores)
      i = 0
      res = foreach (i=1:num_cores,.combine = c) %dopar% {
        library(perms)
        St = Sdiv
        if(i<=rest){
          St= St+1
        }
        .Call(C_get_log_perms_bioassay, as.numeric(t(X[starts[i]:stops[i],])), as.numeric(levels), as.integer(successes), as.integer(trials), as.integer(n), as.integer(num_trials),as.integer(St), as.integer(debug))
      }
      
    }
    
  }
  
  
  
  return(res)
}