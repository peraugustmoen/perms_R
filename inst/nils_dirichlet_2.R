# dirichlet nils example

library(fastperm)
n = 100

set.seed(10002)

ts = seq(-3,3, length.out = n/10)
ts = rep(ts, each=10)

number_of_succ = c(0,0,2,1,4,6,9,10,10,10)
  
I = rep(0, n)
for (i in 1:10) {
  if(number_of_succ[i]>0){
    start = ((i-1)*10+1 )
    stop = ((i)*10 )
    I[start : (start+number_of_succ[i] -1)] = 1
  }
  
}
I = I>0

#I =  (U <= (1/3 * pnorm((ts+2)/0.7) + 2/3 * pnorm((ts-1)/0.7)))

xleft = ts[I]
xright = ts[!I]

nleft = length(xleft)
nright = length(xright)

a = c(rep(-Inf, nleft), xright)
b = c(xleft, rep(Inf, nright))


T = 50002

alpha_param = 1

X = matrix(NA, nrow = T, ncol = n)
# generate samples from prior
for (t in 1:T) {
  X[t,1] = rnorm(1)
  Us = runif(n-1)
  for (i in 2:n) {
    if(Us[i-1] <= alpha_param / (alpha_param +i-1)){
      X[t,i] = rnorm(1)
    }else{
      if(i==2){
        X[t,i] = X[t,i-1]
      }else{
        X[t,i] = sample(x = X[t,1:(i-1)],1)
      }
      if(X[t,i]==1){
        print("her")
      }
    }
    
  }
}


# nonzeroperms = rep(TRUE,T)
# for (t in 1:T) {
#   X[t,] = sort(X[t,])
#   for (i in 1:n){
#     if(X[t,i] < a[i] || X[t,i] > b[i]){
#       nonzeroperms[t] = FALSE
#       #print(sprintf("perm must be zero! "))
#     }
#   }
#   
# }

#sum(nonzeroperms)
# nils has now generated T marginal samples from the *marginal* Dirichlet preusess

system.time({logperms = get_log_permanents(X,a,b, debug = FALSE, parallel=TRUE)})




log_ess = 2 *TEST_logsumexp(logperms[logperms>-1]) - TEST_logsumexp(2*logperms[logperms>-1])
ess = exp(log_ess)
ess


#computing marginal likelihood

log_factorial = function(n){
  s = 0
  for (i in 1:n) {
    s = s + log(i)
  }
  return(s)
}

ml = TEST_logsumexp(logperms[logperms>-1]) - log_factorial(100)  - log(T)
for (i in 1:10) {
  ml = ml + log(choose(10,number_of_succ[i]))
}






system.time({logperms_seq = get_log_permanents(X,a,b, debug = FALSE, parallel=FALSE)})










sum(abs(logperms - logperms_seq))
## E(X|data)
rowzum = rowSums(X)



