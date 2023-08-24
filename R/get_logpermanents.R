#' @useDynLib perms C_get_log_perms
#' @import doParallel
#' @import parallel
#' @import foreach
#' @title get_log_perms
#' @description Computes log permanents associated with simulated latent variables.
#' Each row of the S x n matrix X contains a random sample of size n from the data model. 
#' If there is only a single covariate, then the
#' observed data are represented as (t,y), where t is the observed
#' values of the covariate and y is the vector of indicator variables.
#' If there are more covariates or the problem is phrased as binary
#' classification (see Section 5 in [1]), then t is an S x n matrix
#' since the threshold values change in each iteration. The function returns
#' a vector of log permanents corresponding to each sample in X.
#' @param X A matrix of dimension S x n, in which each row contains a sample from
#' the data model.
#' @param tt Either: A vector of length n containing the observed values of the covariate,
#' Or: A matrix of dimension S x n (if there are several covariates).
#' @param y A vector of length n indicating whether x_i<=t_i for each i in the observed data.
#' @param n Sample size.
#' @param S Number of samples from the data model. That is, the number of iterations in the estimator.
#' @param debug If \code{TRUE}, debug information is printed. 
#' @param parallel If \code{TRUE}, computation is run on several cores
#' @param num_cores (Optional) Specifies the number of cores to use if \code{parallel = TRUE}
#' @return Numpy array of log permanents,each element associated to the corresponding row in X.
#' A zero valued permanent is indicated by a -1.
#' @examples 
#' library(perms)
#' set.seed(1996)
#' n = 100
#' t = seq(0, 1, length.out=n)
#' y = c(rep(0, n/2), rep(1, n/2))
#' S = 200
#' X = matrix(runif(n*S),nrow = S, ncol = n)
#'
#' logperms = get_log_perms(X, t, y, n, S, debug = FALSE, parallel = FALSE, 
#'                          num_cores = NULL)
#' logML = get_log_ML(logperms, n, S, FALSE)
#' logML
#' @references
#' [1] Christensen, D (2023). Inference for Bayesian nonparametric models with binary response data via permutation counting. Bayesian Analysis, Advance online publication, DOI: 10.1214/22-BA1353.
#' @export
get_log_perms= function(X, tt, y, n, S, debug=FALSE, parallel = TRUE, num_cores = NULL){
  
  # X is n times T
  constant_ts = 1
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
    if(length(tt)!= n){
      stop("t must have length n when S = 1")
    }
    if(length(y)!= n){
      stop("y must have length n")
    }

    
    res = .Call(C_get_log_perms, (X[]), as.numeric(tt[]), as.integer(y[]), as.integer(n), as.integer(S), as.integer(debug), as.integer(constant_ts))
    
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
    if(length(y)!= n){
      stop("y must have length n")
    }
    
    if(is.matrix(tt)){
      constant_ts = 0
      if(dim(tt)[1] != S || dim(tt)[2] != n){
        stop("t must either be a vector of length n or an S x n matrix")
      }
      tt = t(tt)
    }else{
      if(length(tt)!=n){
        stop("t must have length n when given as a vector")
      }
      tt = as.numeric(tt)
    }
    
    
    if(!parallel){
      res = .Call(C_get_log_perms, t(X[]), tt, as.integer(y[]), as.integer(n), as.integer(S), as.integer(debug), as.integer(constant_ts))
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
        .Call(C_get_log_perms, t(X[starts[i]:stops[i],]), tt, as.integer(y[]), as.integer(n), as.integer(St), as.integer(debug),as.integer(constant_ts))
      }
      
    }
   
  }
  else{
    stop("Input X must be an S x n matrix, or a vector of length n if S = 1")
  }
  
 
  return(res)
}