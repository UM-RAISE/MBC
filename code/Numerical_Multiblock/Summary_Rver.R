# Multi-block Parameter Calibration 

# #########################################################################
# Description: Summary 
# Date: 2022.8.27
# #########################################################################

rm(list=ls(all=TRUE))
# 
# setwd("C:/Research/Multi-block Simulation Parameter Calibration")
# rslt_path = "C:/Research/Multi-block Simulation Parameter Calibration/Result"

# #########################################################################
# Example 3: L2 calibration (Table 6)
# #########################################################################

table6 = c()

load(file = paste0("/root/capsule/results/L2_result.RData"))

r1 = res_param_mean[idx,]
r2 = res_param_sd[idx,]
r3 = res_mse_mean[idx,]
r4 = res_mse_sd[idx,]

table6 = c(r1[1],r2[1],r1[2],r2[2],r1[3],r2[3],r3[1],r4[1],r3[2],r4[2],r3[3],r4[3])

print(table6)

# #########################################################################
# Example 4: high-dimensional case (Table 13)
# #########################################################################

table13 = c()

naming = "highdim"
# load(file = paste0(rslt_path,"/M_BC_result_",naming,".RData"))
load(file = paste0("/root/capsule/results/M_BC_result_",naming,".RData"))

idx = which.min(func_train_mean_store[,2])

result = c()
aa = matrix(theta_mean_store[idx,])
bb = matrix(theta_std_store[idx,])
cc = matrix(func_test_mean_store[idx,2:(B+1)])
dd = matrix(func_test_std_store[idx,2:(B+1)])
tmp = rhs_band_list[[idx]]
idx.not.na = apply((is.na(tmp) == 0),1,prod)
tmp1 = tmp[as.logical(idx.not.na),]
ee = matrix(apply(tmp1,2,mean))
ff = matrix(apply(tmp1,2,sd))

LB = theta_mle_list[[idx]][as.logical(idx.not.na),] - rhs_band_list[[idx]][as.logical(idx.not.na),]
UB = theta_mle_list[[idx]][as.logical(idx.not.na),] + rhs_band_list[[idx]][as.logical(idx.not.na),]

coverage = c()
for (i in c(1:B)) {
  coverage[i] = 100*sum( (i >= LB[,i]) & (i <= UB[,i]) )/length(LB[,i])
}
gg = matrix(coverage)

result = cbind(aa,bb,cc,dd,ee,ff,gg)
time_rslt = cbind(time_rslt,time)

table13 = result

print(table13)

# #########################################################################
# Example 5: block dominance case (Table 14)
# #########################################################################

table14 = c()
time_rslt = c()

naming = "dominant"
# load(file = paste0(rslt_path,"/M_BC_result_",naming,".RData"))
load(file = paste0("/root/capsule/results/M_BC_result_",naming,".RData"))

idx = which.min(func_train_mean_store[,2])

result = c()
aa = matrix(theta_mean_store[idx,])
bb = matrix(theta_std_store[idx,])
cc = matrix(func_test_mean_store[idx,2:(B+1)])
dd = matrix(func_test_std_store[idx,2:(B+1)])
tmp = rhs_band_list[[idx]]
idx.not.na = apply((is.na(tmp) == 0),1,prod)
tmp1 = tmp[as.logical(idx.not.na),]
ee = matrix(apply(tmp1,2,mean))
ff = matrix(apply(tmp1,2,sd))

LB = theta_mle_list[[idx]][as.logical(idx.not.na),] - rhs_band_list[[idx]][as.logical(idx.not.na),]
UB = theta_mle_list[[idx]][as.logical(idx.not.na),] + rhs_band_list[[idx]][as.logical(idx.not.na),]

coverage = c()
for (i in c(1:B)) {
  coverage[i] = 100*sum( (i >= LB[,i]) & (i <= UB[,i]) )/length(LB[,i])
}
gg = matrix(coverage)

result = cbind(aa,bb,cc,dd,ee,ff,gg)
time_rslt = cbind(time_rslt,time)

table14 = rbind(table14,result)

#######################################

naming = "dominant"
# load(file = paste0(rslt_path,"/H_BC_result_",naming,".RData"))
load(file = paste0("/root/capsule/results/H_BC_result_",naming,".RData"))

idx = which.min(func_train_mean_store[,2])

result = c()
aa = matrix(theta_mean_store[idx,])
bb = matrix(theta_std_store[idx,])
cc = matrix(func_test_mean_store[idx,2:(B+1)])
dd = matrix(func_test_std_store[idx,2:(B+1)])
tmp = rhs_band_list[[idx]]
idx.not.na = apply((is.na(tmp) == 0),1,prod)
tmp1 = tmp[as.logical(idx.not.na),]
ee = matrix(apply(tmp1,2,mean))
ff = matrix(apply(tmp1,2,sd))

LB = theta_mle_list[[idx]][as.logical(idx.not.na),] - rhs_band_list[[idx]][as.logical(idx.not.na),]
UB = theta_mle_list[[idx]][as.logical(idx.not.na),] + rhs_band_list[[idx]][as.logical(idx.not.na),]

coverage = c()
for (i in c(1:B)) {
  coverage[i] = 100*sum( (i >= LB[,i]) & (i <= UB[,i]) )/length(LB[,i])
}
gg = matrix(coverage)

result = cbind(aa,bb,cc,dd,ee,ff,gg)
time_rslt = cbind(time_rslt,time)

table14 = rbind(table14,result)

#######################################

naming = "dominant"
# load(file = paste0(rslt_path,"/S_BC_result_",naming,".RData"))
load(file = paste0("/root/capsule/results/S_BC_result_",naming,".RData"))

idx = which.min(func_train_mean_store[,2])

result = c()
aa = matrix(theta_mean_store[idx,])
bb = matrix(theta_std_store[idx,])
cc = matrix(func_test_mean_store[idx,2:(B+1)])
dd = matrix(func_test_std_store[idx,2:(B+1)])
tmp = rhs_band_list[[idx]]
idx.not.na = apply((is.na(tmp) == 0),1,prod)
tmp1 = tmp[as.logical(idx.not.na),]
ee = matrix(apply(tmp1,2,mean))
ff = matrix(apply(tmp1,2,sd))

LB = theta_mle_list[[idx]][as.logical(idx.not.na),] - rhs_band_list[[idx]][as.logical(idx.not.na),]
UB = theta_mle_list[[idx]][as.logical(idx.not.na),] + rhs_band_list[[idx]][as.logical(idx.not.na),]

coverage = c()
for (i in c(1:B)) {
  coverage[i] = 100*sum( (i >= LB[,i]) & (i <= UB[,i]) )/length(LB[,i])
}
gg = matrix(coverage)

result = cbind(aa,bb,cc,dd,ee,ff,gg)
time_rslt = cbind(time_rslt,time)

table14 = rbind(table14,result)

print(table14)


