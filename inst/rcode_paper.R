#library(devtools)
#install_github("peraugustmoen/perms_R")

library(perms)
set.seed(1996)
n = 100
t = seq(0, 1, length.out=n)
y = c(rep(0, n/2), rep(1, n/2))
S = 20000
X = matrix(runif(n*S),nrow = S, ncol = n)

logperms = get_log_perms(X, t, y, n, S, debug = FALSE, parallel = FALSE, num_cores = NULL)

num_nonzero_perms = sum(logperms > -1)

logML = get_log_ML(logperms, n, S, FALSE)
logML


# exact answer
answer = 0
for(i in 1:n){
  if (y[i] == 1){
    answer = answer + log(t[i])
  }

  else{
    answer = answer + log(1 - t[i])
  }
}
print(answer)

## ten runs
runs = 10
logml_vec = rep(NA, runs)
times_vec = rep(NA, runs)
times_par_vec = rep(NA, runs)
num_nonzero_perms_vec = rep(NA, runs)


for (i in 1:runs) {
  X = matrix(runif(n*S),nrow = S, ncol = n)
  times_vec[i] = system.time({
    logperms = get_log_perms(X, t, y, n, S, debug = FALSE, parallel = FALSE, num_cores = NULL)
    logml_vec[i] = get_log_ML(logperms, n, S, FALSE)
    num_nonzero_perms_vec[i] = sum(logperms > -1)
    })[3]
  times_par_vec[i] = system.time({
    logperms = get_log_perms(X, t, y, n, S, debug = FALSE, parallel = TRUE, num_cores = NULL)
    #logml_vec[i] = get_log_ML(logperms, n, S, FALSE)
    #num_nonzero_perms_vec[i] = sum(logperms > -1)
  })[3]
}
proportion = num_nonzero_perms_vec / S*100
  
logml_vec
times_vec
times_par_vec
proportion


## Dirichlet toy model

set.seed(1996)
n = 500
num_trials = 10
levels = seq(-1, 1, length.out = num_trials)

trials = rep(n %/% num_trials, num_trials)
successes = c(10, 26, 10, 20, 20, 19, 29, 24, 31, 33)

S = 100000
alpha = 1.0

get_X = function(S,n,alpha,seed){
  set.seed(seed)
  X = matrix(0, nrow = S, ncol = n)
  for (s in 1:S) {
    X[s,1] = rnorm(1)
    for (i in 2:n) {
      u = runif(1)
      if(u < (alpha/(alpha+i-1))){
        X[s,i] = rnorm(1)
      }else{
        if(i==2){
          X[s,i] = X[s,1]
        }else{
          X[s,i] = sample(X[s, 1:(i-1)],size=1)
        }
      }
      
    }
    
  }
  return(X)
}
           

runs = 10
logml_vec = rep(NA, runs)
times_vec = rep(NA, runs)
times_par_vec = rep(NA, runs)
num_nonzero_perms_vec = rep(NA, runs)

for(i in 1:runs){
  seed = 1000+i
  set.seed(seed)
  X = get_X(S, n, alpha, seed)
  times_vec[i] = system.time({
    logperms = get_log_perms_bioassay(X, levels, successes, trials, n, num_trials, S,debug=FALSE,parallel = FALSE)
    logml_vec[i] = get_log_ML_bioassay(logperms, successes, trials, n, num_trials, S)
    num_nonzero_perms_vec[i] = sum(logperms > -1)
  })[3]
  
  times_par_vec[i] = system.time({
    logperms = get_log_perms_bioassay(X, levels, successes, trials, n, num_trials, S,parallel = TRUE, debug=FALSE)
    #logml_vec[i] = get_log_ML_bioassay(logperms, successes, trials, n, num_trials, S)
    #num_nonzero_perms_vec[i] = sum(logperms > -1)
  })[3]
  if(abs(get_log_ML_bioassay(logperms, successes, trials, n, num_trials, S) - logml_vec[i])>1e-10){
    print("ERROR THEY ARE DIFFERENT")
  }
  print(logml_vec[i])
  
}
proportion = num_nonzero_perms_vec / S*100

logml_vec
times_vec
times_par_vec
proportion

round(logml_vec, digits = 4)
round(proportion, digits = 3)
round(times_vec, digits = 2)
