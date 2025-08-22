% Multi-block Parameter Calibration

% #########################################################################
% Description: H-BC Algorithm
% Date: 2021.9.22
% #########################################################################

clear all; close all; clc; 

% #########################################################################
% 1. Problem instance
% - EX1: convex + perfect model (case_no = 1)
% - EX2: nonconvex + perfect model (case_no = 2)
% - EX3: nonconvex + imperfect model (case_no = 3)
% #########################################################################

% problem setting 
case_no = 1; 

n = 1000; % data size 
nTry = 10; %200; % number of experiments
dt = 10^(-8); % finite difference
epsilon = 10^(-4); % convergence tolerance
max_iter = 2000;

I = @(t) (t<pi);

Ia = @(t) (t<(2/5)*pi);
Ib = @(t) (t>=(2/5)*pi).*(t<(4/5)*pi);
Ic = @(t) (t>=(4/5)*pi).*(t<(6/5)*pi);
Id = @(t) (t>=(6/5)*pi).*(t>=(8/5)*pi);
Ie = @(t) (t>=(8/5)*pi);

% true parameter
Ttheta0 = -1; 
Ttheta1 = 10; 
Ttheta2 = 5; 

% starting points using permutation 
start = [-10 -10 -10; -10 -10 20; -10 20 -10; -10 20 20; 20 -10 -10; 20 -10 20; 20 20 -10; 20 20 20]; % EX1
% start = [-10 9 4; -10 9 6; -10 11 4; -10 11 6; 10 9 4; 10 9 6; 10 11 4; 10 11 6]; % EX2,3 

% noise level 
noise_level = 0.01;

% #########################################################################
% 2. Algorithm
% #########################################################################

theta_mean_store = []; theta_std_store = []; % calibrated values
F_mean_store = []; F_std_store = []; % training loss 
F_test_mean_store = []; F_test_std_store = []; % testing loss
F_eval_mean_store = []; F_eval_std_store = []; % function evaluations 
outer_iter_mean_store = []; outer_iter_std_store = []; % outer + inner iterations
theta_mle_cell = {};
sigma2_mle_cell = {};
N_cell = {};
rhs_band_cell = {};
tic;

for j = 1:size(start,1) 

    rng('default')
    
    start_theta_rslt = [];
    start_F_rslt = [];
    start_F_test_rslt = [];
    start_F_eval_rlst = [];
    start_outer_iter_rlst = [];
    start_N_rslt = [];
    start_band_rslt = [];

    for iTry = 1:nTry        
    
    disp(['experiment:' num2str(iTry)]) 
    
    if (case_no == 1) || (case_no == 2)
        % training data ~ Unif(0,2*pi)
        x_init = 2*pi.*rand(n,1);
        x0 = x_init(x_init<(4/5)*pi);
        x1 = x_init(logical((x_init>=(2/5)*pi).*(x_init<(8/5)*pi))); 
        x2 = x_init(x_init>=(6/5)*pi);
        
        noise_init = sqrt(noise_level).*randn(n,1);
        noise0 = noise_init(x_init<(4/5)*pi);
        noise1 = noise_init(logical((x_init>=(2/5)*pi).*(x_init<(8/5)*pi))); 
        noise2 = noise_init(x_init>=(6/5)*pi);
        
        % testing data ~ Unif(0,2*pi)
        x_test_init = 2*pi.*rand(n,1);
        x0_test = x_test_init(x_test_init<(4/5)*pi);
        x1_test = x_test_init(logical((x_test_init>=(2/5)*pi).*(x_test_init<(8/5)*pi))); 
        x2_test = x_test_init(x_test_init>=(6/5)*pi);
        
        noise_test_init = sqrt(noise_level).*randn(n,1);
        noise0_test = noise_test_init(x_test_init<(4/5)*pi);
        noise1_test = noise_test_init(logical((x_test_init>=(2/5)*pi).*(x_test_init<(8/5)*pi)));
        noise2_test = noise_test_init(x_test_init>=(6/5)*pi);
    else       
        % training data ~ Unif(0,2*pi)
        x_init = 2*pi.*rand(n,1);
        x0 = x_init; 
        x1 = x_init(x_init<=pi);
        x2 = x_init(x_init>pi);
        
        noise_init = sqrt(noise_level).*randn(n,1);
        noise0 = noise_init;
        noise1 = noise_init(x_init<pi);
        noise2 = noise_init(x_init>=pi);
        
        % testing data ~ Unif(0,2*pi)
        x_test_init = 2*pi.*rand(n,1);
        x0_test = x_test_init;
        x1_test = x_test_init(x_test_init<pi);
        x2_test = x_test_init(x_test_init>=pi);
        
        noise_test_init = sqrt(noise_level).*randn(n,1);
        noise0_test = noise_test_init;
        noise1_test = noise_test_init(x_test_init<pi);
        noise2_test = noise_test_init(x_test_init>=pi);
    end 
    
    % problem  
    if case_no == 1
        fpt = @(x,theta0,theta1,theta2,noiset) Ia(x).*exp(x./10).*sin(x) + Ib(x).*exp(x./10).*cos(x) + Ic(x).*exp(x./10).*cos(x) + Id(x).*exp(x./10).*cos(x) + Ie(x).*exp(x./10).*sin(x) + noiset;
        fp0 = @(x,theta0,theta1,theta2,noise0) Ia(x).*exp(x./10).*sin(x) + Ib(x).*exp(x./10).*cos(x) + noise0; 
        fp1 = @(x,theta0,theta1,theta2,noise1)                             Ib(x).*exp(x./10).*cos(x) + Ic(x).*exp(x./10).*cos(x) + Id(x).*exp(x./10).*cos(x) + noise1;
        fp2 = @(x,theta0,theta1,theta2,noise2)                                                                                     Id(x).*exp(x./10).*cos(x) + Ie(x).*exp(x./10).*sin(x) + noise2; 

        fct = @(x,theta0,theta1,theta2) Ia(x).*exp(x./10).*sin(x) + Ib(x).*exp(x./10).*cos(x) + Ic(x).*exp(x./10).*cos(x) + Id(x).*exp(x./10).*cos(x) + Ie(x).*exp(x./10).*sin(x)   - Ia(x).*(1/2)*(theta0+1).*sqrt(x.^2-x+1) - Ib(x).*(1/10).*(theta0+1)*(theta1-10).*sqrt(x.^2-x+1) - Ic(x).*(1/10)*(theta1-10).*sqrt(x.^2-x+1) - Id(x).*(1/10)*(theta1-10)*(theta2-5).*sqrt(x.^2-x+1) - Ie(x).*(1/20).*(theta2-5).*sqrt(x.^2-x+1);
        fc0 = @(x,theta0,theta1,theta2) Ia(x).*exp(x./10).*sin(x) + Ib(x).*exp(x./10).*cos(x)                                                                                       - Ia(x).*(1/2)*(theta0+1).*sqrt(x.^2-x+1) - Ib(x).*(1/10).*(theta0+1)*(theta1-10).*sqrt(x.^2-x+1); 
        fc1 = @(x,theta0,theta1,theta2)                             Ib(x).*exp(x./10).*cos(x) + Ic(x).*exp(x./10).*cos(x) + Id(x).*exp(x./10).*cos(x)                                                                         - Ib(x).*(1/10).*(theta0+1)*(theta1-10).*sqrt(x.^2-x+1) - Ic(x).*(1/10)*(theta1-10).*sqrt(x.^2-x+1) - Id(x).*(1/10)*(theta1-10)*(theta2-5).*sqrt(x.^2-x+1);
        fc2 = @(x,theta0,theta1,theta2)                                                                                     Id(x).*exp(x./10).*cos(x) + Ie(x).*exp(x./10).*sin(x)                                                                                                                                                 - Id(x).*(1/10).*(theta1-10)*(theta2-5).*sqrt(x.^2-x+1) - Ie(x).*(1/20).*(theta2-5).*sqrt(x.^2-x+1);            

    elseif case_no == 2
        fpt = @(x,theta0,theta1,theta2,noiset) Ia(x).*exp(x./10).*sin(x) + Ib(x).*exp(x./10).*cos(x) + Ic(x).*exp(x./10).*cos(x) + Id(x).*exp(x./10).*cos(x) + Ie(x).*exp(x./10).*sin(x) + noiset;
        fp0 = @(x,theta0,theta1,theta2,noise0) Ia(x).*exp(x./10).*sin(x) + Ib(x).*exp(x./10).*cos(x) + noise0; 
        fp1 = @(x,theta0,theta1,theta2,noise1)                             Ib(x).*exp(x./10).*cos(x) + Ic(x).*exp(x./10).*cos(x) + Id(x).*exp(x./10).*cos(x) + noise1;
        fp2 = @(x,theta0,theta1,theta2,noise2)                                                                                     Id(x).*exp(x./10).*cos(x) + Ie(x).*exp(x./10).*sin(x) + noise2; 
            
        fct = @(x,theta0,theta1,theta2) Ia(x).*exp(x./10).*sin(x) + Ib(x).*exp(x./10).*cos(x) + Ic(x).*exp(x./10).*cos(x) + Id(x).*exp(x./10).*cos(x) + Ie(x).*exp(x./10).*sin(x)   - Ia(x).*(2)*(theta0+1).*(cos(theta0.*x./10)) - Ib(x).*(1/2).*(theta0+1)*(theta1-10).*(cos(theta1.*x./20)) - Ic(x).*(1/2)*(theta1-10).*(cos(theta1.*x./10)) - Id(x).*(1/2).*(theta1-10)*(theta2-5).*(cos(theta1.*x./20)) - Ie(x).*(1/2).*(theta2-5).*(cos(theta2.*x./10));
        fc0 = @(x,theta0,theta1,theta2) Ia(x).*exp(x./10).*sin(x) + Ib(x).*exp(x./10).*cos(x)                                                                                       - Ia(x).*(2)*(theta0+1).*(cos(theta0.*x./10)) - Ib(x).*(1/2).*(theta0+1)*(theta1-10).*(cos(theta1.*x./20)); 
        fc1 = @(x,theta0,theta1,theta2)                             Ib(x).*exp(x./10).*cos(x) + Ic(x).*exp(x./10).*cos(x) + Id(x).*exp(x./10).*cos(x)                                                                             - Ib(x).*(1/2).*(theta0+1)*(theta1-10).*(cos(theta1.*x./20)) - Ic(x).*(1/2)*(theta1-10).*(cos(theta1.*x./10)) - Id(x).*(1/2).*(theta1-10)*(theta2-5).*(cos(theta1.*x./20));
        fc2 = @(x,theta0,theta1,theta2)                                                                                     Id(x).*exp(x./10).*cos(x) + Ie(x).*exp(x./10).*sin(x)                                                                                                                                                               - Id(x).*(1/2).*(theta1-10)*(theta2-5).*(cos(theta1.*x./20)) - Ie(x).*(1/2).*(theta2-5).*(cos(theta2.*x./10));
       
    elseif case_no == 3
        fpt = @(x,theta0,theta1,theta2,noiset) exp(x./10) + I(x).*exp(x./10).*sin(x) + (1-I(x)).*exp(x./10).*cos(x) + noiset;
        fp0 = @(x,theta0,theta1,theta2,noise0) exp(x./10) + I(x).*exp(x./10).*sin(x) + (1-I(x)).*exp(x./10).*cos(x) + noise0; 
        fp1 = @(x,theta0,theta1,theta2,noise1) exp(x./10) + exp(x./10).*sin(x) + noise1;
        fp2 = @(x,theta0,theta1,theta2,noise2) exp(x./10) + exp(x./10).*cos(x) + noise2; 
        
        fct = @(x,theta0,theta1,theta2) exp(x./10) + I(x).*exp(x./10).*sin(x) + (1-I(x)).*exp(x./10).*cos(x) - 5.*sqrt(theta0.^2-theta0+1).*(sin(theta0.*x./10)+cos(theta0.*x./10)) - I(x).*(1/2).*sqrt(theta1.^2-theta1+1).*(sin(theta1.*x./5)+cos(theta1.*x./5)) - (1/2).*(1-I(x)).*sqrt(theta2.^2-theta2+1).*(sin(theta2.*x./10)+cos(theta2.*x./10));
        fc0 = @(x,theta0,theta1,theta2) exp(x./10) + I(x).*exp(x./10).*sin(x) + (1-I(x)).*exp(x./10).*cos(x) - 5.*sqrt(theta0.^2-theta0+1).*(sin(theta0.*x./10)+cos(theta0.*x./10)) - I(x).*(1/2).*sqrt(theta1.^2-theta1+1).*(sin(theta1.*x./5)+cos(theta1.*x./5)) - (1/2).*(1-I(x)).*sqrt(theta2.^2-theta2+1).*(sin(theta2.*x./10)+cos(theta2.*x./10));
        fc1 = @(x,theta0,theta1,theta2) exp(x./10) + exp(x./10).*sin(x) - 5.*sqrt(theta0.^2-theta0+1).*(sin(theta0.*x./10)+cos(theta0.*x./10)) - (1/2).*sqrt(theta1.^2-theta1+1).*(sin(theta1.*x./5)+cos(theta1.*x./5)); 
        fc2 = @(x,theta0,theta1,theta2) exp(x./10) + exp(x./10).*cos(x) - 5.*sqrt(theta0.^2-theta0+1).*(sin(theta0.*x./10)+cos(theta0.*x./10)) - (1/2).*sqrt(theta2.^2-theta2+1).*(sin(theta2.*x./10)+cos(theta2.*x./10));    
    end
    
    % objective function
    Ft = @(x,theta0,theta1,theta2,noiset) (1/length(x))*sum((fct(x,theta0,theta1,theta2) - fpt(x,Ttheta0,Ttheta1,Ttheta2,noiset)).^2);
    F0 = @(x,theta0,theta1,theta2,noise0) (1/length(x))*sum((fc0(x,theta0,theta1,theta2) - fp0(x,Ttheta0,Ttheta1,Ttheta2,noise0)).^2);
    F1 = @(x,theta0,theta1,theta2,noise1) (1/length(x))*sum((fc1(x,theta0,theta1,theta2) - fp1(x,Ttheta0,Ttheta1,Ttheta2,noise1)).^2);
    F2 = @(x,theta0,theta1,theta2,noise2) (1/length(x))*sum((fc2(x,theta0,theta1,theta2) - fp2(x,Ttheta0,Ttheta1,Ttheta2,noise2)).^2);

    % gradient
%     gF0 = @(x,theta0,theta1,theta2) (F0(x,theta0+dt,theta1,theta2,noise0) - F0(x,theta0-dt,theta1,theta2,noise0))/(2*dt);
%     gF1 = @(x,theta0,theta1) (F1(x,theta0,theta1+dt,noise1) - F1(x,theta0,theta1-dt,noise1))/(2*dt);
%     gF2 = @(x,theta0,theta2) (F2(x,theta0,theta2+dt,noise2) - F2(x,theta0,theta2-dt,noise2))/(2*dt);      

    d1 = @(x,theta0,theta1,theta2,noiset) (Ft(x,theta0+dt,theta1,theta2,noiset) - Ft(x,theta0-dt,theta1,theta2,noiset))/(2*dt); 
    d2 = @(x,theta0,theta1,theta2,noiset) (Ft(x,theta0,theta1+dt,theta2,noiset) - Ft(x,theta0,theta1-dt,theta2,noiset))/(2*dt);
    d3 = @(x,theta0,theta1,theta2,noiset) (Ft(x,theta0,theta1,theta2+dt,noiset) - Ft(x,theta0,theta1,theta2-dt,noiset))/(2*dt);
    
    % for Fisher information matrix
    like0 = @(x,theta0,theta1,theta2,noiset,mse) -(1/2)*log(2*pi) -(1/2)*log(mse) -(1/(2*mse))*(fct(x,theta0,theta1,theta2) - fpt(x,Ttheta0,Ttheta1,Ttheta2,noiset)).^2;
    
    I00 = @(x,theta0,theta1,theta2,noiset,mse) -(1/length(x))*sum((-like0(x,theta0+2*dt,theta1,theta2,noiset,mse) + 16*like0(x,theta0+dt,theta1,theta2,noiset,mse) - 30*like0(x,theta0,theta1,theta2,noiset,mse) + 16*like0(x,theta0-dt,theta1,theta2,noiset,mse) - like0(x,theta0-2*dt,theta1,theta2,noiset,mse))/(12*dt^2));
    I01 = @(x,theta0,theta1,theta2,noiset,mse) -(1/length(x))*sum((like0(x,theta0+dt,theta1+dt,theta2,noiset,mse) - like0(x,theta0+dt,theta1-dt,theta2,noiset,mse) - like0(x,theta0-dt,theta1+dt,theta2,noiset,mse) + like0(x,theta0-dt,theta1-dt,theta2,noiset,mse))/(4*dt^2) );
    I02 = @(x,theta0,theta1,theta2,noiset,mse) -(1/length(x))*sum((like0(x,theta0+dt,theta1,theta2+dt,noiset,mse) - like0(x,theta0+dt,theta1,theta2-dt,noiset,mse) - like0(x,theta0-dt,theta1,theta2+dt,noiset,mse) + like0(x,theta0-dt,theta1,theta2-dt,noiset,mse))/(4*dt^2) );
   
    I10 = @(x,theta0,theta1,theta2,noiset,mse) -(1/length(x))*sum((like0(x,theta0+dt,theta1+dt,theta2,noiset,mse) - like0(x,theta0-dt,theta1+dt,theta2,noiset,mse) - like0(x,theta0+dt,theta1-dt,theta2,noiset,mse) + like0(x,theta0-dt,theta1-dt,theta2,noiset,mse))/(4*dt^2) );
    I11 = @(x,theta0,theta1,theta2,noiset,mse) -(1/length(x))*sum((-like0(x,theta0,theta1+2*dt,theta2,noiset,mse) + 16*like0(x,theta0,theta1+dt,theta2,noiset,mse) - 30*like0(x,theta0,theta1,theta2,noiset,mse) + 16*like0(x,theta0,theta1-dt,theta2,noiset,mse) - like0(x,theta0,theta1-2*dt,theta2,noiset,mse))/(12*dt^2) );
    I12 = @(x,theta0,theta1,theta2,noiset,mse) -(1/length(x))*sum((like0(x,theta0,theta1+dt,theta2+dt,noiset,mse) - like0(x,theta0,theta1+dt,theta2-dt,noiset,mse) - like0(x,theta0,theta1-dt,theta2+dt,noiset,mse) + like0(x,theta0,theta1-dt,theta2-dt,noiset,mse))/(4*dt^2) );
    
    I20 = @(x,theta0,theta1,theta2,noiset,mse) -(1/length(x))*sum((like0(x,theta0+dt,theta1,theta2+dt,noiset,mse) - like0(x,theta0-dt,theta1,theta2+dt,noiset,mse) - like0(x,theta0+dt,theta1,theta2-dt,noiset,mse) + like0(x,theta0-dt,theta1,theta2-dt,noiset,mse))/(4*dt^2) );
    I21 = @(x,theta0,theta1,theta2,noiset,mse) -(1/length(x))*sum((like0(x,theta0,theta1+dt,theta2+dt,noiset,mse) - like0(x,theta0,theta1-dt,theta2+dt,noiset,mse) - like0(x,theta0,theta1+dt,theta2-dt,noiset,mse) + like0(x,theta0,theta1-dt,theta2-dt,noiset,mse))/(4*dt^2) );
    I22 = @(x,theta0,theta1,theta2,noiset,mse) -(1/length(x))*sum((-like0(x,theta0,theta1,theta2+2*dt,noiset,mse) + 16*like0(x,theta0,theta1,theta2+dt,noiset,mse) - 30*like0(x,theta0,theta1,theta2,noiset,mse) + 16*like0(x,theta0,theta1,theta2-dt,noiset,mse) - like0(x,theta0,theta1,theta2-2*dt,noiset,mse))/(12*dt^2) );  
        
    % initialization
    Itheta = [start(j,1) start(j,2) start(j,3)]';
    theta = Itheta;

    k = 1; % number of iteration
      
    f_eval_F = 0; 
    
    while 1       
        [tau0,f_eval_F] = back_with(theta,x_init,noise_init,case_no,f_eval_F);
        gF_new = [d1(x_init,theta(1),theta(2),theta(3),noise_init);d2(x_init,theta(1),theta(2),theta(3),noise_init);d3(x_init,theta(1),theta(2),theta(3),noise_init)];
        theta = theta - tau0*gF_new;     
        f_eval_F = f_eval_F + 6;
        
        % training loss
        Fta = Ft(x_init,theta(1),theta(2),theta(3),noise_init);
        F0a = F0(x0,theta(1),theta(2),theta(3),noise0); f_eval_F = f_eval_F + 1;
        F1a = F1(x1,theta(1),theta(2),theta(3),noise1);
        F2a = F2(x2,theta(1),theta(2),theta(3),noise2);        
        
        % testing loss
        Ftb = Ft(x_test_init,theta(1),theta(2),theta(3),noise_test_init);
        F0b = F0(x0_test,theta(1),theta(2),theta(3),noise0_test); f_eval_F = f_eval_F + 1;
        F1b = F1(x1_test,theta(1),theta(2),theta(3),noise1_test);
        F2b = F2(x2_test,theta(1),theta(2),theta(3),noise2_test);
        
        % convergence criteria    
        new_F0 = F0a; new_theta = theta; 
        
        if k ~= 1
%             delta = abs(new_F0-old_F0)/abs(old_F0);
            delta = norm(gF_new,2); f_eval_F = f_eval_F + 6;
            if (delta<epsilon) || (k>max_iter) 
                break;
            end
        end

        old_F0 = new_F0; old_theta = new_theta;
        k = k + 1;
    end  
    tempF = [F0a F1a F2a Fta]; 
    tempF_test = [F0b F1b F2b Ftb];  
    tempF_N = [length(x0) length(x1) length(x2)];
    Imat = [I00(x_init,theta(1),theta(2),theta(3),noise_init,Fta) I01(x_init,theta(1),theta(2),theta(3),noise_init,Fta) I02(x_init,theta(1),theta(2),theta(3),noise_init,Fta);
            I10(x_init,theta(1),theta(2),theta(3),noise_init,Fta) I11(x_init,theta(1),theta(2),theta(3),noise_init,Fta) I12(x_init,theta(1),theta(2),theta(3),noise_init,Fta);
            I20(x_init,theta(1),theta(2),theta(3),noise_init,Fta) I21(x_init,theta(1),theta(2),theta(3),noise_init,Fta) I22(x_init,theta(1),theta(2),theta(3),noise_init,Fta)];
    invI = inv(Imat); 
    temp_band = [1.96*sqrt(invI(1,1))/sqrt(length(x_init)) ...
                 1.96*sqrt(invI(2,2))/sqrt(length(x_init)) ...
                 1.96*sqrt(invI(3,3))/sqrt(length(x_init))] 
    start_theta_rslt = [start_theta_rslt; theta'];
    start_F_rslt = [start_F_rslt; tempF];
    start_F_test_rslt = [start_F_test_rslt; tempF_test];
    start_F_eval_rlst = [start_F_eval_rlst; f_eval_F];
    start_outer_iter_rlst = [start_outer_iter_rlst; k];
    start_N_rslt = [start_N_rslt; tempF_N];
    start_band_rslt = [start_band_rslt; temp_band];
    
    disp(['= result: ' num2str(theta(1)) ', ' num2str(theta(2)) ', ' num2str(theta(3)) ]) 
    
    end
%     [~,idx] = min(start_F_rslt(:,1));
%     best_theta_rslt = start_theta_rslt(:,idx);
%     best_F_rslt = start_F_rslt(idx,:);
%     
%     theta_rslt = [theta_rslt; best_theta_rslt'];
%     F_rslt = [F_rslt; best_F_rslt];

    theta_mean_store = [theta_mean_store; mean(start_theta_rslt,1)];
    theta_std_store = [theta_std_store; std(start_theta_rslt,1)];
    F_mean_store = [F_mean_store; mean(start_F_rslt,1)];
    F_std_store = [F_std_store; std(start_F_rslt,1)];
    F_test_mean_store = [F_test_mean_store; mean(start_F_test_rslt,1)];
    F_test_std_store = [F_test_std_store; std(start_F_test_rslt,1)];
    F_eval_mean_store = [F_eval_mean_store; mean(start_F_eval_rlst,1)];
    F_eval_std_store = [F_eval_std_store; std(start_F_eval_rlst,1)];
    outer_iter_mean_store = [outer_iter_mean_store; mean(start_outer_iter_rlst,1)];
    outer_iter_std_store = [outer_iter_std_store; std(start_outer_iter_rlst,1)];

    theta_mle_cell{j,1} = start_theta_rslt;
    sigma2_mle_cell{j,1} = start_F_rslt;
    N_cell{j,1} = start_N_rslt;
    rhs_band_cell{j,1} = start_band_rslt;
    
end       

theta_mean_store
theta_std_store
F_mean_store
F_std_store
F_test_mean_store
F_test_std_store

% F_eval_mean_store
% F_eval_std_store % function evaluations 
% outer_iter_mean_store
% outer_iter_std_store % outer interations
time = toc

% select the result whose traning loss is minimum 
% [~,idx] = min(F_mean_store,[],1);
% idx = idx(1);
% result = [theta_mean_store(idx,:) theta_std_store(idx,:) ...
%             F_test_mean_store(idx,:) F_test_std_store(idx,:) ...
%             F_eval_mean_store(idx,:) 0 0 F_eval_std_store(idx,:) 0 0 ...
%             outer_iter_mean_store(idx,:) outer_iter_std_store(idx,:) ...
%             outer_iter_mean_store(idx,:) outer_iter_std_store(idx,:) time ...
%             real(mean(rhs_band_cell{idx,1},1))];
%  
% real(mean(rhs_band_cell{idx,1},1))
% 
% LB = theta_mle_cell{idx,1} - real(rhs_band_cell{idx,1});
% UB = theta_mle_cell{idx,1} + real(rhs_band_cell{idx,1});
% ci0 = 100*sum( (Ttheta0 >= LB(:,1)) & (Ttheta0 <= UB(:,1)) )/nTry;
% ci1 = 100*sum( (Ttheta1 >= LB(:,2)) & (Ttheta1 <= UB(:,2)) )/nTry;
% ci2 = 100*sum( (Ttheta2 >= LB(:,3)) & (Ttheta2 <= UB(:,3)) )/nTry;
% coverage = [ci0 ci1 ci2]

% coverage_rslt = [];
% half_confidence_length = [];
% zero_idx = find(prod(real(rhs_band_cell{idx,1}),2) == 0);
% 
% LB = theta_mle_cell{idx,1} - rhs_band_cell{idx,1};
% UB = theta_mle_cell{idx,1} + rhs_band_cell{idx,1};
% 
% rhs_band_cell{idx,1}(zero_idx,:) = [];
% LB(zero_idx,:) = [];
% UB(zero_idx,:) = [];
% half_confidence_length = [half_confidence_length; mean(rhs_band_cell{idx,1},1)];
% 
% ci0 = 100*sum( (Ttheta0 >= LB(:,1)) & (Ttheta0 <= UB(:,1)) )/length(LB);
% ci1 = 100*sum( (Ttheta1 >= LB(:,2)) & (Ttheta1 <= UB(:,2)) )/length(LB);
% ci2 = 100*sum( (Ttheta2 >= LB(:,3)) & (Ttheta2 <= UB(:,3)) )/length(LB);
% coverage = [ci0 ci1 ci2]
    
% #########################################################################
% 3. Save file
% #########################################################################

save(['/results/EX' num2str(case_no) '_H_BC.mat'])

