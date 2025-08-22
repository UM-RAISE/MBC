# Multi-block Parameter Calibration 

# #########################################################################
# Description: H-BC Algorithm
# Date: 2022.8.27
# #########################################################################

rm(list=ls(all=TRUE))

# setwd("C:/Research/Multi-block Simulation Parameter Calibration")
# rslt_path = "C:/Research/Multi-block Simulation Parameter Calibration/Result"

library(lhs)

set.seed(1)

# #########################################################################
# 1. Problem instance
# #########################################################################

naming = "dominant"

# algorithmic parameters 
n = 50 # data size 
n_test = 100
B = 3 # number of block 
omega = 0 # shrinkage factor (D=0: overlap, D=1/2: no overlap)
nTry = 5 #200 # number of experiments
h = 10^(-8) # finite difference
epsilon = 10^(-4) # convergence tolerance
max_iter = 200 # 2000 
max_x = 4*pi
D = 6 # number of subblock 

# backtracking linesearch parameters
alpha_bar = 1
rho = 0.5
c1 = 10^(-4)

# starting point: 0~(B+1) lhs 
# init = randomLHS(8,B)*(B+1)
init = (randomLHS(8,B)-0.5)*20

noise_level = 0.01

# true process
zeta <- function(x,B,D) {
  z = ( x < 2*max_x/(2*D+1) )*exp(x/10)*sin(x) +
    ( (x >= max_x/(2*D+1)) & (x < 2*D*max_x/(2*D+1)) )*exp(x/10)*cos(x) +
    ( x >= 2*max_x*(2*D-1)/(2*D+1) )*exp(x/10)*sin(x)    
  return(z) 
}

# computer model 
eta <- function(x,theta,B,D) {
  e1 = (1/5)*(theta[1] - 1)*sqrt(x^2-x+1)
  e2 = (theta[1] - 1)*(theta[2] - 2)*sqrt(x^2-x+1)
  e3 = (1/5)*(theta[2] - 2)*sqrt(x^2-x+1)
  e4 = (theta[2] - 2)*(theta[3] - 3)*sqrt(x^2-x+1) 
  e5 = (1/5)*(theta[3] - 3)*sqrt(x^2-x+1)
  e = e1 + e2 + e3 + e4 + e5 + zeta(x,B,D)  
  return(e)
}

# true process by block
zetaB <- function(x,b,B,D) {
  zb = 0
  if (b/2==1) {
    zb = exp(x/10)*cos(x) 
  } else {
    zb = exp(x/10)*sin(x) 
  }
  return(zb)
}

# computer model by block
etaB <- function(x,theta,b,B,D) {
  if (b==1) {
    eb = (1/5)*(theta[1] - 1)*sqrt(x^2-x+1) + (theta[1] - 1)*(theta[2] - 2)*sqrt(x^2-x+1)
  } else if (b==B) {
    eb = (theta[2] - 2)*(theta[3] - 3)*sqrt(x^2-x+1) + (1/5)*(theta[3] - 3)*sqrt(x^2-x+1)
  } else {
    eb = (theta[1] - 1)*(theta[2] - 2)*sqrt(x^2-x+1) + (1/5)*(theta[2] - 2)*sqrt(x^2-x+1) + (theta[2] - 2)*(theta[3] - 3)*sqrt(x^2-x+1)
  }
  eb = eb + zetaB(x,b,B,D)
  return(eb)
}

# objective function 
func <- function(x,theta,B,D,noise) {
  f = (1/length(x))*sum( ( (zeta(x,B,D)+noise) - eta(x,theta,B,D) )^2 )
  return(f)
}

# objective function by block 
funcB <- function(x,theta,b,B,D,noise) {
  f = (1/length(x))*sum( ( (zetaB(x,b,B,D)+noise) - etaB(x,theta,b,B,D) )^2 )
  return(f)
}

# gradient 
grad <- function(x,theta,B,D,noise) {
  g = matrix(0,B,1)
  for (i in c(1:B)) {
    g[i,1] = (func(x,theta+h*tabulate(i,nbins=B),B,D,noise) - func(x,theta-h*tabulate(i,nbins=B),B,D,noise))/(2*h)
  }
  return(g)
}

# gradient by block
gradB <- function(x,theta,b,B,D,noise) {
  g = numeric(0)
  g = (funcB(x,theta+h*tabulate(b,nbins=B),b,B,D,noise) - funcB(x,theta-h*tabulate(b,nbins=B),b,B,D,noise))/(2*h)
  return(g)
}

# Hessian 
like = function(x,theta,B,D,noise,mse) {
  l = -(1/2)*log(2*pi) -(1/2)*log(mse) -(1/(2*mse))*( (zeta(x,B,D)+noise) - eta(x,theta,B,D) )^2
  return(l)
}
Hlike = function(x,theta,B,D,noise,mse) {
  H = matrix(0,B,B)
  for (i in c(1:B)) {
    for (j in c(1:B)) {
      if (i == j) {
        H[i,j] = -(1/length(x))*sum( (-like(x,theta+2*h*tabulate(i,nbins=B),B,D,noise,mse) + 16*like(x,theta+h*tabulate(i,nbins=B),B,D,noise,mse) -
                                        30*like(x,theta,B,D,noise,mse) + 16*like(x,theta-h*tabulate(i,nbins=B),B,D,noise,mse) -
                                        like(x,theta-2*h*tabulate(i,nbins=B),B,D,noise,mse))/(12*h^2) )
      } else {
        H[i,j] = -(1/length(x))*sum( (like(x,theta+h*tabulate(i,nbins=B)+h*tabulate(j,nbins=B),B,D,noise,mse) - 
                                        like(x,theta+h*tabulate(i,nbins=B)-h*tabulate(j,nbins=B),B,D,noise,mse) - 
                                        like(x,theta-h*tabulate(i,nbins=B)+h*tabulate(j,nbins=B),B,D,noise,mse) + 
                                        like(x,theta-h*tabulate(i,nbins=B)-h*tabulate(j,nbins=B),B,D,noise,mse))/(4*h^2) )
      }
    }
  }
  H = H + diag(B)*10^(-8)
  return(H)
} 

# Hessian by block
likeB = function(x,theta,b,B,D,noise,mse) {
  l = -(1/2)*log(2*pi) -(1/2)*log(mse) -(1/(2*mse))*( (zetaB(x,b,B,D)+noise) - etaB(x,theta,b,B,D) )^2
  return(l)
}
HlikeB = function(x,theta,b,B,D,noise,mse) {
  H = -(1/length(x))*sum((-likeB(x,theta+2*h*tabulate(b,nbins=B),b,B,D,noise,mse) + 16*likeB(x,theta+h*tabulate(b,nbins=B),b,B,D,noise,mse) -
                            30*likeB(x,theta,b,B,D,noise,mse) + 16*likeB(x,theta-h*tabulate(b,nbins=B),b,B,D,noise,mse) -
                            likeB(x,theta-2*h*tabulate(b,nbins=B),b,B,D,noise,mse))/(12*h^2))
  return(H)
} 

# #########################################################################
# 2. Algorithm
# #########################################################################

theta_mean_store = matrix(0,nrow(init),B)
theta_std_store = matrix(0,nrow(init),B)
func_train_mean_store = matrix(0,nrow(init),B+1)
func_train_std_store = matrix(0,nrow(init),B+1)
func_test_mean_store = matrix(0,nrow(init),B+1)
func_test_std_store = matrix(0,nrow(init),B+1)

theta_mle_list = list()
sigma2_mle_list = list()
n_list = list()
rhs_band_list = list()

tic <- proc.time()[3]
for (j in c(1:nrow(init))) {
  
  x_train_list = list()
  x_test_list = list() 
  noise_train_list = list()
  noise_test_list = list() 
  
  init_theta_rslt = matrix(0,nTry,B)
  init_func_train_rslt = matrix(0,nTry,B+1)
  init_func_test_rslt = matrix(0,nTry,B+1)
  init_n_rslt = matrix(0,nTry,B)
  init_band_rslt = matrix(0,nTry,B) 
  
  for (iTry in c(1:nTry)) {
    
    cat("experiment:", iTry, "\n")  
    theta = init[j,] # initialization
    
    set.seed(iTry)  
    
    x_train = runif(n,min=0,max=max_x)
    
    temp = x_train*( x_train < 2*max_x/(2*D+1) )
    x_train_list[[1]] = temp[temp!=0] 
    temp = x_train*( (x_train >= max_x/(2*D+1)) & (x_train < (2*D)*max_x/(2*D+1)) )
    x_train_list[[2]] = temp[temp!=0] 
    temp = x_train*( x_train >= (2*D-1)*max_x/(2*D+1) )
    x_train_list[[3]] = temp[temp!=0] 
    
    noise_train = rnorm(n,mean=0,sd=sqrt(noise_level))
    
    temp = noise_train*( x_train < 2*max_x/(2*D+1) )
    noise_train_list[[1]] = temp[temp!=0] 
    temp = noise_train*( (x_train >= max_x/(2*D+1)) & (x_train < (2*D)*max_x/(2*D+1)) )
    noise_train_list[[2]] = temp[temp!=0] 
    temp = noise_train*( x_train >= (2*D-1)*max_x/(2*D+1) )
    noise_train_list[[3]] = temp[temp!=0]     
    
    set.seed(iTry+2022)
    x_test = runif(n_test,min=0,max=max_x)
    
    temp = x_test*( x_test < 2*max_x/(2*D+1) )
    x_test_list[[1]] = temp[temp!=0] 
    temp = x_test*( (x_test >= max_x/(2*D+1)) & (x_test < (2*D)*max_x/(2*D+1)) )
    x_test_list[[2]] = temp[temp!=0]     
    temp = x_test*( x_test >= (2*D-1)*max_x/(2*D+1) )
    x_test_list[[3]] = temp[temp!=0]       
    
    noise_test = rnorm(n_test,mean=0,sd=sqrt(noise_level))
    
    temp = noise_test*( x_test < 2*max_x/(2*D+1) )
    noise_test_list[[1]] = temp[temp!=0] 
    temp = noise_test*( (x_test >= max_x/(2*D+1)) & (x_test < (2*D)*max_x/(2*D+1)) )
    noise_test_list[[2]] = temp[temp!=0] 
    temp = noise_test*( x_test >= (2*D-1)*max_x/(2*D+1) )
    noise_test_list[[3]] = temp[temp!=0] 
    
    k = 1
    
    while (TRUE) {
      
      # linesearch 
      alpha = alpha_bar 
      g = grad(x_train,theta,B,D,noise_train)
      while ( func(x_train,theta-alpha*g,B,D,noise_train) > func(x_train,theta,B,D,noise_train) - c1*alpha*norm(g,type="2")^2 ) {
        alpha = alpha*rho
      }
      # GD update 
      theta = theta-alpha*g
      
      # training and testing loss
      func_tr = numeric(0)
      func_tr[1] = func(x_train,theta,B,D,noise_train)
      func_te = numeric(0)
      func_te[1] = func(x_test,theta,B,D,noise_test)
      n_assign = numeric(0)
      
      for (t in c(1:B)) {
        n_assign[t] = length(x_train_list[[t]])
        x_tr = x_train_list[[t]]
        noise_tr = noise_train_list[[t]]
        func_tr[t+1] = funcB(x_tr,theta,t,B,D,noise_tr)
        x_te = x_test_list[[t]]
        noise_te = noise_test_list[[t]]
        func_te[t+1] = funcB(x_te,theta,t,B,D,noise_te)
      }
      
      band_length = numeric(0)
      mse = func_tr[2]
      H = Hlike(x_train,theta,B,D,noise_train,mse)  
      invH = solve(H)
      
      for (t in c(1:B)) {
        x_tr = x_train_list[[t]]
        noise_tr = noise_train_list[[t]]
        band_length[t] = 1.96*sqrt(invH[t,t])/sqrt(length(x_train)) 
      }
      
      if (k!=1) {
        delta = norm(g,"2")
        if (delta<epsilon || k>max_iter) {
          break
        }
      }
      
      k = k + 1
    } # iteration 
    
    init_theta_rslt[iTry,] = theta
    init_func_train_rslt[iTry,] = func_tr
    init_func_test_rslt[iTry,] = func_te
    init_n_rslt[iTry,] = n_assign
    init_band_rslt[iTry,] = band_length # right half bandwidth
    
    cat("= result:", theta, "\n")
  } # experiment 
  
  theta_mean_store[j,] = apply(init_theta_rslt,2,mean)
  theta_std_store[j,] = apply(init_theta_rslt,2,sd)
  func_train_mean_store[j,] = apply(init_func_train_rslt,2,mean)
  func_train_std_store[j,] = apply(init_func_train_rslt,2,sd)
  func_test_mean_store[j,] = apply(init_func_test_rslt,2,mean)
  func_test_std_store[j,] = apply(init_func_test_rslt,2,sd)
  
  theta_mle_list[[j]] = init_theta_rslt
  sigma2_mle_list[[j]] = init_func_train_rslt
  n_list[[j]] = init_n_rslt
  rhs_band_list[[j]] = init_band_rslt
  
} # starting point 
toc <- proc.time()[3]
time = toc - tic 
print(time)

# #########################################################################
# 3. Save file
# #########################################################################

# save.image(file = paste0(rslt_path,"/S_BC_result_",naming,".RData"))
# load(file = paste0(rslt_path,"/S_BC_result_",naming,".RData"))

save.image(file = paste0("/root/capsule/results/S_BC_result_",naming,".RData"))
# load(file = paste0("/root/capsule/results/S_BC_result_",naming,".RData"))
