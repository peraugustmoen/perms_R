# dirichlet nils example

library(fastperm)
n = 100

set.seed(10002)

ts = seq(-1,1, length.out = n)

U = runif(n)

I =  (U <= (1/3 * pnorm((ts+2)/0.7) + 2/3 * pnorm((ts-1)/0.7)))

xleft = ts[I]
xright = ts[!I]

nleft = length(xleft)
nright = length(xright)

a = c(rep(-Inf, nleft), xright)
b = c(xleft, rep(Inf, nright))


T = 200

alpha_param = 100

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
  X[t,] = sort(X[t,])
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

logperms



log_ess = 2 *TEST_logsumexp(logperms[logperms>-1]) - TEST_logsumexp(2*logperms[logperms>-1])
ess = exp(log_ess)
ess




system.time({logperms_seq = get_log_permanents(X,a,b, debug = FALSE, parallel=FALSE)})










sum(abs(logperms - logperms_seq))
## E(X|data)
rowzum = rowSums(X)



