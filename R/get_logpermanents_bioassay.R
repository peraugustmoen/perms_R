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
#'  log permanents corresponding to each sample.
#' @param X A matrix of dimension S x n, in which each row contains a sample from
#' the data model.
#' @param levels A vector of length n containing the levels at which trials were conducted.
#' @param successes A vector of length n containing the number of successful trials at each level.
#' @param trials A vector of length n containing the number of trials at each level.
#' @param n Sample size.
#' @param num_trials Number of different trials. 
#' @param S Number of samples from the data model. That is, the number of iterations in the estimator.
#' @param debug If \code{TRUE}, debug information is printed. 
#' @param parallel If \code{TRUE}, computation is run on several cores
#' @param num_cores (Optional) Specifies the number of cores to use if \code{parallel = TRUE}
#' @return Numpy array of log permanents, each element associated to the corresponding row in X.
#' A zero valued permanent is indicated by a -1.
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
#' logperms = get_log_perms_bioassay(X, levels, successes, trials, n, num_trials, S,
#'            debug=FALSE,parallel = FALSE)
#' logml = get_log_ML_bioassay(logperms, successes, trials, n, num_trials, S)
#' 
#' proportion = sum(logperms>-1) / S*100
#' 
#' proportion 
#' logml
#' @references
#' [1] Christensen, D (2023). Inference for Bayesian nonparametric models with binary response data via permutation counting. Bayesian Analysis, Advance online publication, DOI: 10.1214/22-BA1353.
#' @export
get_log_perms_bioassay= function(X, levels, successes, trials, n, num_trials, S, debug=FALSE, parallel = TRUE, num_cores = NULL){
  
  # X is n times T
  if(!is.vector(levels) || !is.vector(successes) || !is.vector(trials)){
    stop("levels, successes, trials must all be vectors")
  }
  if(length(levels)!= num_trials){
    stop("levels must have length num_trials")
  }
  if(length(trials)!= num_trials){
    stop("trials must have length num_trials")
  }
  if(length(successes)!= num_trials){
    stop("successes must have length num_trials")
  }
  if(sum(trials)!= n){
    stop("sum of all elements of trials must be equal to n")
  }
    
  if(is.vector(X)){
    
    if(debug){
      print("X is input as vector")
    }
    if(n != length(X)){
      stop("X must have length n when given as vector")
    }
    if(S != 1){
      stop("When S>1, X must be an S x n matrix")
    }
    
    
    
    
    res = .Call(C_get_log_perms_bioassay, (X[]), as.numeric(levels), as.integer(successes), as.integer(trials), as.integer(n), as.integer(num_trials),as.integer(S), as.integer(debug))
    
  }else if(is.matrix(X)){
    if(debug){
      print("X is input as matrix:")
    }
    
    if(S !=dim(X)[1] || n != dim(X)[2]){
      stop("X must be S x n whenever S>1")
    }
    if(S <= 1){
      stop("When X is a matrix, S must be >1")
    }
    
    
    
    if(!parallel){
      res = .Call(C_get_log_perms_bioassay, t(X[]),as.numeric(levels), as.integer(successes), as.integer(trials), as.integer(n), as.integer(num_trials),as.integer(S), as.integer(debug))
    }else{
      #print("par")
      if(is.null(num_cores)){
        num_cores <- detectCores()
        if(debug){
          print(sprintf("Registered %d cores\n", as.integer(num_cores)))
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
        .Call(C_get_log_perms_bioassay, t(X[starts[i]:stops[i],]), as.numeric(levels), as.integer(successes), as.integer(trials), as.integer(n), as.integer(num_trials),as.integer(St), as.integer(debug))
      }
      
    }
    
  }
  else{
    stop("Input X must be an S x n matrix, or a vector of length n if S = 1")
  }
  
  
  return(res)
}