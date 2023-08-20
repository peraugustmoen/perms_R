library(perms)

for (i in 1:100) {
  n = 3
  
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
  
  
  S = 2000
  
  alpha_param = 100
  
  X = matrix(NA, nrow = S, ncol = n)
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
  
  t2 = matrix(NA, nrow = S, ncol = n)
  
  for(i in 1:S){
    t2[i,] = t[]
  }
  
  
  
  res2 = get_log_permanents(X,t2,y,n,S,FALSE)
  print(res2)
  
  sum(abs(res - res2))
  
  
}

x = c(-3, -0.5, 0.07)
y = c(0,0,1)
t = c(-1,0,1)

get_log_permanents(x,t,y,3,1,FALSE)


x = c(-0.5, 0.5, 2)
y = c(0,1,1)
t = c(-1,0,1)

get_log_permanents(x,t,y,3,1,FALSE)




