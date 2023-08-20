library(perms)
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


S = 200

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
  
}


y = rep(0,n)
t = rep(0.0,n)

for (i in 1:n){
  if (a[i]<=-10000){
    y[i] = 1
    t[i] = b[i]
  }
  else{
    t[i] = a[i]
  }
  
}


res = get_log_permanents(X,t,y,n,S,FALSE)
print(res)

res2 = get_log_permanents(X,t,y,n,S,FALSE,TRUE)
print(res2)

sum(abs(res-res2))

res3 = get_log_permanents(X,t,y,n,S,FALSE,TRUE,6)
print(res3)

sum(abs(res-res3))
## samme resultat som python!

## tester med feil inputs

res = get_log_permanents(X[1:10,],t,y,n,S,FALSE)
res = get_log_permanents(X[,1:10],t,y,n,S,FALSE)
res = get_log_permanents(X,t[1:10],y,n,S,FALSE)
res = get_log_permanents(X,t,y[1:10],n,S,FALSE)
