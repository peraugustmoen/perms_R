#' @useDynLib perms C_wrapper_log_sum_exp

#' @title log_sum_exp
#' @description Computes the log sum exp of a vector. 
#' Given input array = [x_1, ..., x_n], returns x_* + log(exp(x_1 - x_*) + ... + exp(x_n - x_*)),
#' where x_* = max(x_1, ... x_n). Ignores entries
#' with NA value.
#' @param x A vector length n containing the computed log permanents, where a zero permanent
#'  is indicated by a NA value.
#' @param x Input vector. 
#' @return The log-sum-exp of the entries of the input vector. 
#' @examples
#' library(perms)
#' x = c(1,2,3,-1,-1,1)
#' log_sum_exp(x)
#' @export
log_sum_exp= function(x){
  res = .Call(C_wrapper_log_sum_exp, x, as.integer(length(x)))
  return(res)
}
