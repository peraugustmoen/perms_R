library(perms)


N = 10
hh = 1
n = N*hh

set.seed(10002)

levels = seq(-1,1, length.out = N)
ts = rep(levels, each = hh)

U = runif(n)

I =  (U <= (1/3 * pnorm((ts+2)/0.7) + 2/3 * pnorm((ts-1)/0.7)))

xleft = ts[I]
xright = ts[!I]

nleft = length(xleft)
nright = length(xright)

a = c(rep(-Inf, nleft), xright)
b = c(xleft, rep(Inf, nright))

trials = rep(hh, N)
successes = rep(NA, N)
for(i in 1:N){
  successes[i] = sum(I[(hh*(i-1)+1): (hh*i)])
}
S = 100

alpha_param = 100

X = matrix(NA, nrow=S, ncol =n)


# generate samples from prior
for (t in 1:S) {
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

if(S==1){
  X = as.vector(X)
}





y = rep(1, n)
t = rep(1.0, n)

for(i in 1:n){
  if (a[i]<=-10000){
    y[i] = 1
    t[i] = b[i]
  }else{
    t[i] = a[i]
    y[i] = 0
  }
}
  
#print(t)
#print(y)

# standard:
res = get_log_permanents(X,t,y,n,S,FALSE)
#print(res)

# bioassay:
#print(successes)
resb = get_log_permanents_bioassay(X,levels,successes,trials,n,N,S,FALSE)
#print(resb)

print("difference")
print(sum(abs(res-resb)))


logml1 = get_log_ML(res, n,S)
print("ml1:")
print(logml1)
logml2 =get_log_ML_bioassay(resb,successes,trials,n,N,S)
print("ml2:")
print(logml2)




