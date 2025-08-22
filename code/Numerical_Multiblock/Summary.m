% Multi-block Parameter Calibration 

% #########################################################################
% Description: Summary
% Date: 2021.9.22
% #########################################################################

clear all; close all; clc; 

% #########################################################################
% Table 3
% #########################################################################

summary = [];

load(['/results/EX1_M_BC.mat'])
summary = [summary; result];
load(['/results/EX1_H_BC.mat'])
summary = [summary; result];
load(['/results/EX1_S_BC.mat'])
summary = [summary; result];

load(['/results/EX2_M_BC.mat'])
summary = [summary; result];
load(['/results/EX2_H_BC.mat'])
summary = [summary; result];
load(['/results/EX2_S_BC.mat'])
summary = [summary; result];

load(['/results/EX3_M_BC.mat'])
summary = [summary; result];
load(['/results/EX3_H_BC.mat'])
summary = [summary; result];
load(['/results/EX3_S_BC.mat'])
summary = [summary; result];

csvwrite('/results/TableA.txt',M)

% #########################################################################
% Table 7
% #########################################################################

summary = [];

load(['/results/EX1_M_BC.mat'])
summary = [summary; result];
load(['/results/EX1_M_BC_V.mat'])
summary = [summary; result];
load(['/results/EX1_M_BC_VS.mat'])
summary = [summary; result];

load(['/results/EX2_M_BC.mat'])
summary = [summary; result];
load(['/results/EX2_M_BC_V.mat'])
summary = [summary; result];
load(['/results/EX2_M_BC_VS.mat'])
summary = [summary; result];

load(['/results/EX3_M_BC.mat'])
summary = [summary; result];
load(['/results/EX3_M_BC_V.mat'])
summary = [summary; result];
load(['/results/EX3_M_BC_VS.mat'])
summary = [summary; result];

csvwrite('/results/TableB.txt',M)

% #########################################################################
% Table 5 
% #########################################################################

coverage_rslt = [];
half_confidence_length = [];
half_confidence_length_std = [];

load(['/results/EX1_M_BC.mat'])
zero_idx = find(prod(real(rhs_band_cell{idx,1}),2) == 0);

LB = theta_mle_cell{idx,1} - rhs_band_cell{idx,1};
UB = theta_mle_cell{idx,1} + rhs_band_cell{idx,1};

rhs_band_cell{idx,1}(zero_idx,:) = [];
LB(zero_idx,:) = [];
UB(zero_idx,:) = [];
half_confidence_length = [half_confidence_length; mean(rhs_band_cell{idx,1},1)];
half_confidence_length_std = [half_confidence_length_std; std(rhs_band_cell{idx,1},1)];

ci0 = 100*sum( (Ttheta0 >= LB(:,1)) & (Ttheta0 <= UB(:,1)) )/length(LB);
ci1 = 100*sum( (Ttheta1 >= LB(:,2)) & (Ttheta1 <= UB(:,2)) )/length(LB);
ci2 = 100*sum( (Ttheta2 >= LB(:,3)) & (Ttheta2 <= UB(:,3)) )/length(LB);
coverage = [ci0 ci1 ci2];
coverage_rslt = [coverage_rslt; coverage];

load(['/results/EX1_H_BC.mat'])
zero_idx = find(prod(real(rhs_band_cell{idx,1}),2) == 0);

LB = theta_mle_cell{idx,1} - rhs_band_cell{idx,1};
UB = theta_mle_cell{idx,1} + rhs_band_cell{idx,1};

rhs_band_cell{idx,1}(zero_idx,:) = [];
LB(zero_idx,:) = [];
UB(zero_idx,:) = [];
half_confidence_length = [half_confidence_length; mean(rhs_band_cell{idx,1},1)];
half_confidence_length_std = [half_confidence_length_std; std(rhs_band_cell{idx,1},1)];

ci0 = 100*sum( (Ttheta0 >= LB(:,1)) & (Ttheta0 <= UB(:,1)) )/length(LB);
ci1 = 100*sum( (Ttheta1 >= LB(:,2)) & (Ttheta1 <= UB(:,2)) )/length(LB);
ci2 = 100*sum( (Ttheta2 >= LB(:,3)) & (Ttheta2 <= UB(:,3)) )/length(LB);
coverage = [ci0 ci1 ci2];
coverage_rslt = [coverage_rslt; coverage];

load(['/results/EX1_S_BC.mat'])
zero_idx = find(prod(real(rhs_band_cell{idx,1}),2) == 0);

LB = theta_mle_cell{idx,1} - rhs_band_cell{idx,1};
UB = theta_mle_cell{idx,1} + rhs_band_cell{idx,1};

rhs_band_cell{idx,1}(zero_idx,:) = [];
LB(zero_idx,:) = [];
UB(zero_idx,:) = [];
half_confidence_length = [half_confidence_length; mean(rhs_band_cell{idx,1},1)];
half_confidence_length_std = [half_confidence_length_std; std(rhs_band_cell{idx,1},1)];

ci0 = 100*sum( (Ttheta0 >= LB(:,1)) & (Ttheta0 <= UB(:,1)) )/length(LB);
ci1 = 100*sum( (Ttheta1 >= LB(:,2)) & (Ttheta1 <= UB(:,2)) )/length(LB);
ci2 = 100*sum( (Ttheta2 >= LB(:,3)) & (Ttheta2 <= UB(:,3)) )/length(LB);
coverage = [ci0 ci1 ci2];
coverage_rslt = [coverage_rslt; coverage];

load(['/results/EX2_M_BC.mat'])
zero_idx = find(prod(real(rhs_band_cell{idx,1}),2) == 0);

LB = theta_mle_cell{idx,1} - rhs_band_cell{idx,1};
UB = theta_mle_cell{idx,1} + rhs_band_cell{idx,1};

rhs_band_cell{idx,1}(zero_idx,:) = [];
LB(zero_idx,:) = [];
UB(zero_idx,:) = [];
half_confidence_length = [half_confidence_length; mean(rhs_band_cell{idx,1},1)];
half_confidence_length_std = [half_confidence_length_std; std(rhs_band_cell{idx,1},1)];

ci0 = 100*sum( (Ttheta0 >= LB(:,1)) & (Ttheta0 <= UB(:,1)) )/length(LB);
ci1 = 100*sum( (Ttheta1 >= LB(:,2)) & (Ttheta1 <= UB(:,2)) )/length(LB);
ci2 = 100*sum( (Ttheta2 >= LB(:,3)) & (Ttheta2 <= UB(:,3)) )/length(LB);
coverage = [ci0 ci1 ci2];
coverage_rslt = [coverage_rslt; coverage];

load(['/results/EX2_H_BC.mat'])
zero_idx = find(prod(real(rhs_band_cell{idx,1}),2) == 0);

LB = theta_mle_cell{idx,1} - rhs_band_cell{idx,1};
UB = theta_mle_cell{idx,1} + rhs_band_cell{idx,1};

rhs_band_cell{idx,1}(zero_idx,:) = [];
LB(zero_idx,:) = [];
UB(zero_idx,:) = [];
half_confidence_length = [half_confidence_length; mean(rhs_band_cell{idx,1},1)];
half_confidence_length_std = [half_confidence_length_std; std(rhs_band_cell{idx,1},1)];

ci0 = 100*sum( (Ttheta0 >= LB(:,1)) & (Ttheta0 <= UB(:,1)) )/length(LB);
ci1 = 100*sum( (Ttheta1 >= LB(:,2)) & (Ttheta1 <= UB(:,2)) )/length(LB);
ci2 = 100*sum( (Ttheta2 >= LB(:,3)) & (Ttheta2 <= UB(:,3)) )/length(LB);
coverage = [ci0 ci1 ci2];
coverage_rslt = [coverage_rslt; coverage];

load(['/results/EX2_S_BC.mat'])
zero_idx = find(prod(real(rhs_band_cell{idx,1}),2) == 0);

LB = theta_mle_cell{idx,1} - rhs_band_cell{idx,1};
UB = theta_mle_cell{idx,1} + rhs_band_cell{idx,1};

rhs_band_cell{idx,1}(zero_idx,:) = [];
LB(zero_idx,:) = [];
UB(zero_idx,:) = [];
half_confidence_length = [half_confidence_length; mean(rhs_band_cell{idx,1},1)];
half_confidence_length_std = [half_confidence_length_std; std(rhs_band_cell{idx,1},1)];

ci0 = 100*sum( (Ttheta0 >= LB(:,1)) & (Ttheta0 <= UB(:,1)) )/length(LB);
ci1 = 100*sum( (Ttheta1 >= LB(:,2)) & (Ttheta1 <= UB(:,2)) )/length(LB);
ci2 = 100*sum( (Ttheta2 >= LB(:,3)) & (Ttheta2 <= UB(:,3)) )/length(LB);
coverage = [ci0 ci1 ci2];
coverage_rslt = [coverage_rslt; coverage];

load(['/results/EX3_M_BC.mat'])
zero_idx = find(prod(real(rhs_band_cell{idx,1}),2) == 0);

LB = theta_mle_cell{idx,1} - rhs_band_cell{idx,1};
UB = theta_mle_cell{idx,1} + rhs_band_cell{idx,1};

rhs_band_cell{idx,1}(zero_idx,:) = [];
LB(zero_idx,:) = [];
UB(zero_idx,:) = [];
half_confidence_length = [half_confidence_length; mean(rhs_band_cell{idx,1},1)];
half_confidence_length_std = [half_confidence_length_std; std(rhs_band_cell{idx,1},1)];

ci0 = 100*sum( (Ttheta0 >= LB(:,1)) & (Ttheta0 <= UB(:,1)) )/length(LB);
ci1 = 100*sum( (Ttheta1 >= LB(:,2)) & (Ttheta1 <= UB(:,2)) )/length(LB);
ci2 = 100*sum( (Ttheta2 >= LB(:,3)) & (Ttheta2 <= UB(:,3)) )/length(LB);
coverage = [ci0 ci1 ci2];
coverage_rslt = [coverage_rslt; coverage];

load(['/results/EX3_H_BC.mat'])
zero_idx = find(prod(real(rhs_band_cell{idx,1}),2) == 0);

LB = theta_mle_cell{idx,1} - rhs_band_cell{idx,1};
UB = theta_mle_cell{idx,1} + rhs_band_cell{idx,1};

rhs_band_cell{idx,1}(zero_idx,:) = [];
LB(zero_idx,:) = [];
UB(zero_idx,:) = [];
half_confidence_length = [half_confidence_length; mean(rhs_band_cell{idx,1},1)];
half_confidence_length_std = [half_confidence_length_std; std(rhs_band_cell{idx,1},1)];

ci0 = 100*sum( (Ttheta0 >= LB(:,1)) & (Ttheta0 <= UB(:,1)) )/length(LB);
ci1 = 100*sum( (Ttheta1 >= LB(:,2)) & (Ttheta1 <= UB(:,2)) )/length(LB);
ci2 = 100*sum( (Ttheta2 >= LB(:,3)) & (Ttheta2 <= UB(:,3)) )/length(LB);
coverage = [ci0 ci1 ci2];
coverage_rslt = [coverage_rslt; coverage];

load(['/results/EX3_S_BC.mat'])
zero_idx = find(prod(real(rhs_band_cell{idx,1}),2) == 0);

LB = theta_mle_cell{idx,1} - rhs_band_cell{idx,1};
UB = theta_mle_cell{idx,1} + rhs_band_cell{idx,1};

rhs_band_cell{idx,1}(zero_idx,:) = [];
LB(zero_idx,:) = [];
UB(zero_idx,:) = [];
half_confidence_length = [half_confidence_length; mean(rhs_band_cell{idx,1},1)];
half_confidence_length_std = [half_confidence_length_std; std(rhs_band_cell{idx,1},1)];

ci0 = 100*sum( (Ttheta0 >= LB(:,1)) & (Ttheta0 <= UB(:,1)) )/length(LB);
ci1 = 100*sum( (Ttheta1 >= LB(:,2)) & (Ttheta1 <= UB(:,2)) )/length(LB);
ci2 = 100*sum( (Ttheta2 >= LB(:,3)) & (Ttheta2 <= UB(:,3)) )/length(LB);
coverage = [ci0 ci1 ci2];
coverage_rslt = [coverage_rslt; coverage];

csvwrite('/results/TableC.txt',M)