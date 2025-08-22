# Multi-block Parameter Calibration 

# #########################################################################
# Description: M-BC Algorithm
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

naming = "highdim"

# algorithmic parameters 
n = 50 # training data size 
n_test = 100 # test data size
B = 3 # number of block 
omega = 0 # shrinkage factor (omega=-1/2: overlap, omega=1/2: no overlap)
nTry = 5 # 200 # number of experiments
h = 10^(-8) # finite difference
epsilon = 10^(-4) # convergence tolerance
max_iter_out = 200 # 2000 
max_iter_in = 100 # 1000 
max_x = 2*pi

# backtracking linesearch parameters
alpha_bar = 1
rho = 0.5
c1 = 10^(-4)

# starting point: 0~(B+1) lhs 
# init = randomLHS(8,B)*(B+1)
init = (randomLHS(8,B)-0.5)*10

noise_level = 0.01

# true process
zeta <- function(x,B,omega) {
  z1 = 0; z2 = 0
  for (t in c(1:floor( (B+1)/2)) ) {
    z1 = z1 + ( (x >= ( max_x*(4*t-4)/(2*B+1) + max_x*omega/(2*B+1) )) & (x < ( max_x*(4*t-1)/(2*B+1) - max_x*omega/(2*B+1) )) )*exp(x/10)*sin(x) 
  }
  for (t in c(1:floor(B/2)) ) {
    z2 = z2 + ( (x >= ( max_x*(4*t-6)/(2*B+1) + max_x*omega/(2*B+1) )) & (x < ( max_x*(4*t-3)/(2*B+1) - max_x*omega/(2*B+1) )) )*exp(x/10)*cos(x) 
  }
  z = z1 + z2 
  return(z) 
}

# computer model 
eta <- function(x,theta,B,omega) {
  e1 = 0; e2 = 0; e3 = 0; e4 = 0
  for (t in c(1:B)) {
    e1 = e1 + (1/1)*(theta[t] - t)*sqrt(x^2-x+1)*( (x >= max_x*(2*t-1)/(2*B+1)) & (x < max_x*(2*t)/(2*B+1)) )
  }
  for (t in c(1:(B-1))) {
    e2 = e2 + (1/1)*(theta[t] - t)*(theta[t+1] - (t+1))*sqrt(x^2-x+1)*( (x >= ( max_x*(2*t)/(2*B+1) + max_x*omega/(2*B+1) )) & (x < ( max_x*(2*t+1)/(2*B+1) - max_x*omega/(2*B+1) )) )  
  }
  e3 = (1/1)*(theta[1] - 1)*sqrt(x^2-x+1)*( (x >= max_x*omega/(2*B+1)) & (x < max_x/(2*B+1))  )
  e4 = (1/1)*(theta[B] - B)*sqrt(x^2-x+1)*( (x >= max_x - max_x*omega/(2*B+1)) & (x < max_x) )
  e = e1 + e2 + e3 + e4 + zeta(x,B,omega)
  return(e)
}

# true process by block
zetaB <- function(x,b,B,omega) {
  zb = 0
  if (b/2==1) {
    zb = exp(x/10)*cos(x) 
  } else {
    zb = exp(x/10)*sin(x) 
  }
  return(zb)
}

# computer model by block
etaB <- function(x,theta,b,B,omega) {
  eb = 0
  if (b==1) {
    eb = (1/1)*(theta[1] - 1)*sqrt(x^2-x+1)*( (x >= max_x*omega/(2*B+1)) & (x < 2*max_x/(2*B+1))  ) + 
         (1/1)*(theta[1] - 1)*(theta[2] - 2)*sqrt(x^2-x+1)*( (x >= 2*max_x/(2*B+1)) & (x < (3-omega)*max_x/(2*B+1))  )
  } else if (b==B) {
    eb = (1/1)*(theta[B-1] - (B-1))*(theta[B] - B)*sqrt(x^2-x+1)*( (x >= (2*B-2+omega)*max_x/(2*B+1)) & (x < (2*B-1)*max_x/(2*B+1))  ) + 
         (1/1)*(theta[B] - B)*sqrt(x^2-x+1)*( (x >= (2*B-1)*max_x/(2*B+1)) & (x < (2*B+1-omega)*max_x/(2*B+1))  )
  } else {
    eb = (1/1)*(theta[b-1] - (b-1))*(theta[b] - b)*sqrt(x^2-x+1)*( (x >= (2*b-2+omega)*max_x/(2*B+1)) & (x < (2*b-1)*max_x/(2*B+1))  ) + 
         (1/1)*(theta[b] - b)*sqrt(x^2-x+1)*( (x >= (2*b-1)*max_x/(2*B+1)) & (x < (2*b)*max_x/(2*B+1))  ) + 
         (1/1)*(theta[b] - b)*(theta[b+1] - (b+1))*sqrt(x^2-x+1)*( (x >= (2*b)*max_x/(2*B+1)) & (x < (2*b+1-omega)*max_x/(2*B+1))  )
  }
  eb = eb + zetaB(x,b,B,omega)
  return(eb)
}

# objective function 
func <- function(x,theta,B,omega,noise) {
  f = (1/length(x))*sum( ( (zeta(x,B,omega)+noise) - eta(x,theta,B,omega) )^2 )
  return(f)
}

# objective function by block 
funcB <- function(x,theta,b,B,omega,noise) {
  f = (1/length(x))*sum( ( (zetaB(x,b,B,omega)+noise) - etaB(x,theta,b,B,omega) )^2 )
  return(f)
}

# gradient 
grad <- function(x,theta,B,omega,noise) {
  g = matrix(0,B,1)
  for (i in c(1:B)) {
    g[i,1] = (func(x,theta+h*tabulate(i,nbins=B),B,omega,noise) - func(x,theta-h*tabulate(i,nbins=B),B,omega,noise))/(2*h)
  }
  return(g)
}

# gradient by block
gradB <- function(x,theta,b,B,omega,noise) {
  g = numeric(0)
  g = (funcB(x,theta+h*tabulate(b,nbins=B),b,B,omega,noise) - funcB(x,theta-h*tabulate(b,nbins=B),b,B,omega,noise))/(2*h)
  return(g)
}

# Hessian 
like = function(x,theta,B,omega,noise,mse) {
  l = -(1/2)*log(2*pi) -(1/2)*log(mse) -(1/(2*mse))*( (zeta(x,B,omega)+noise) - eta(x,theta,B,omega) )^2
  return(l)
}
Hlike = function(x,theta,B,omega,noise,mse) {
  H = numeric(0)
  for (i in c(1:B)) {
    H[i] = -(1/length(x))*sum((-like(x,theta+2*h*tabulate(i,nbins=B),B,omega,noise,mse) + 16*like(x,theta+h*tabulate(i,nbins=B),B,omega,noise,mse) -
                                 30*like(x,theta,B,omega,noise,mse) + 16*like(x,theta-h*tabulate(i,nbins=B),B,omega,noise,mse) -
                                 like(x,theta-2*h*tabulate(i,nbins=B),B,omega,noise,mse))/(12*h^2))
  }
  return(H)
} 

# Hessian by block
likeB = function(x,theta,b,B,omega,noise,mse) {
  l = -(1/2)*log(2*pi) -(1/2)*log(mse) -(1/(2*mse))*( (zetaB(x,b,B,omega)+noise) - etaB(x,theta,b,B,omega) )^2
  return(l)
}
HlikeB = function(x,theta,b,B,omega,noise,mse) {
  H = -(1/length(x))*sum((-likeB(x,theta+2*h*tabulate(b,nbins=B),b,B,omega,noise,mse) + 16*likeB(x,theta+h*tabulate(b,nbins=B),b,B,omega,noise,mse) -
                            30*likeB(x,theta,b,B,omega,noise,mse) + 16*likeB(x,theta-h*tabulate(b,nbins=B),b,B,omega,noise,mse) -
                            likeB(x,theta-2*h*tabulate(b,nbins=B),b,B,omega,noise,mse))/(12*h^2))
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
for (j in c(1:8)) { # nrow(init)
  
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
    for (t in c(1:B)) {
      temp = x_train*( (x_train >= max_x*(2*t-2)/(2*B+1) + max_x*omega/(2*B+1)) & (x_train < max_x*(2*t+1)/(2*B+1) - max_x*omega/(2*B+1)) )
      x_train_list[[t]] = temp[temp!=0] 
    }
    noise_train = rnorm(n,mean=0,sd=sqrt(noise_level))
    for (t in c(1:B)) {
      temp = noise_train*( (x_train >= max_x*(2*t-2)/(2*B+1) + max_x*omega/(2*B+1)) & (x_train < max_x*(2*t+1)/(2*B+1) - max_x*omega/(2*B+1)) )
      noise_train_list[[t]] = temp[temp!=0] 
    }
    
    set.seed(iTry+2022)
    x_test = runif(n_test,min=0,max=max_x)
    for (t in c(1:B)) {
      temp = x_test*( (x_test >= max_x*(2*t-2)/(2*B+1) + max_x*omega/(2*B+1)) & (x_test < max_x*(2*t+1)/(2*B+1) - max_x*omega/(2*B+1)) )
      x_test_list[[t]] = temp[temp!=0] 
    }
    noise_test = rnorm(n_test,mean=0,sd=sqrt(noise_level))
    for (t in c(1:B)) {
      temp = noise_test*( (x_test >= max_x*(2*t-2)/(2*B+1) + max_x*omega/(2*B+1)) & (x_test < max_x*(2*t+1)/(2*B+1) - max_x*omega/(2*B+1)) )
      noise_test_list[[t]] = temp[temp!=0] 
    }
    
    k = 1
    
    init_func_train_rslt_temp = matrix(0,1,B+1)    
    while (TRUE) {
      
      n_assign = numeric(0)
      func_tr_temp = numeric(0)
      inner_index = numeric(0)

      for (t in c(1:B)) {
        
        # cat("block:",t,"\n")
        
        m = 1
        n_assign[t] = length(x_train_list[[t]])
        
        while (TRUE) {
          
          # linesearch 
          alpha = alpha_bar 
          x_tr = x_train_list[[t]]
          noise_tr = noise_train_list[[t]]
          g = gradB(x_tr,theta,t,B,omega,noise_tr)
          while ( funcB(x_tr,theta-alpha*g*tabulate(t,nbins=B),t,B,omega,noise_tr) > funcB(x_tr,theta,t,B,omega,noise_tr) - c1*alpha*norm(g,type="2")^2 ) {
            alpha = alpha*rho
          }
          # GD update 
          theta = theta-alpha*g*tabulate(t,nbins=B)
          
          if (m != 1) {
            delta = norm(g,"2")
            if (delta<epsilon || m>max_iter_in) {
              break
            }
          }
          m = m + 1
          ##### store F for function trajectory  
          func_tr_temp[1] = func(x_train,theta,B,omega,noise_train)
          x_tr_temp = x_train_list[[t]]
          noise_tr_temp = noise_train_list[[t]]
          func_tr_temp[t+1] = funcB(x_tr_temp,theta,t,B,omega,noise_tr_temp)
          init_func_train_rslt_temp = rbind(init_func_train_rslt_temp,func_tr_temp)
          #####         
        }
        inner_index = cbind(inner_index,m)
      } # inner iteration 

      # training loss
      func_tr = numeric(0)
      func_tr[1] = func(x_train,theta,B,omega,noise_train)
      for (t in c(1:B)) {
        x_tr = x_train_list[[t]]
        noise_tr = noise_train_list[[t]]
        func_tr[t+1] = funcB(x_tr,theta,t,B,omega,noise_tr)
      }
      
      # testing loss
      func_te = numeric(0)
      func_te[1] = func(x_test,theta,B,omega,noise_test)
      for (t in c(1:B)) {
        x_te = x_test_list[[t]]
        noise_te = noise_test_list[[t]]
        func_te[t+1] = funcB(x_te,theta,t,B,omega,noise_te)
      }
      
      band_length = numeric(0)
      mse = func_tr[2]
      for (t in c(1:B)) {
        x_tr = x_train_list[[t]]
        noise_tr = noise_train_list[[t]]
        band_length[t] = 1.96*sqrt(HlikeB(x_tr,theta,t,B,omega,noise_tr,mse)^(-1))/sqrt(length(x_tr)) 
      }
      
      if (k != 1) {
        xi = 0
        for (t in c(1:B)){
          xi = max(xi, norm(gradB(x_train_list[[t]],theta,t,B,omega,noise_train[[t]]),"2") )
        }
        if (xi<epsilon || k>max_iter_out) {
          break
        }
      }
      
      k = k + 1
    } # outer iteration 
    
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

# save.image(file = paste0(rslt_path,"/M_BC_result_",naming,".RData"))
# load(file = paste0(rslt_path,"/M_BC_result_",naming,".RData"))

save.image(file = paste0("/root/capsule/results/M_BC_result_",naming,".RData"))
# load(file = paste0("/root/capsule/results/M_BC_result_",naming,".RData"))