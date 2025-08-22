% Multi-block Parameter Calibration

% #########################################################################
% Description: M-BC-VS
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
case_no = 2; 

n = 1000; % data size 
nTry = 10; %200; % number of experiments
dt = 10^(-8); % finite difference
epsilon = 10^(-4); % convergence tolerance
max_iter_out = 2000; 
max_iter_in = 1000; 

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
% start = [-10 -10 -10; -10 -10 20; -10 20 -10; -10 20 20; 20 -10 -10; 20 -10 20; 20 20 -10; 20 20 20]; % EX1
start = [-10 9 4; -10 9 6; -10 11 4; -10 11 6; 10 9 4; 10 9 6; 10 11 4; 10 11 6]; % EX2,3

% noise level 
noise_level = 0.01;

% #########################################################################
% 2. Algorithm
% #########################################################################

theta_mean_store = []; theta_std_store = []; % calibrated values
F_mean_store = []; F_std_store = []; % training loss 
F_test_mean_store = []; F_test_std_store = []; % testing loss
F_eval_mean_store = []; F_eval_std_store = []; % function evaluations 
total_iter_mean_store = []; total_iter_std_store = []; % outer + inner iterations
outer_iter_mean_store = []; outer_iter_std_store = []; % outer interations
tic;

for j = 1:size(start,1) 

    rng('default')
    
    start_theta_rslt = [];
    start_F_rslt = [];
    start_F_test_rslt = [];
    start_F_eval_rlst = [];
    start_total_iter_rlst = [];
    start_outer_iter_rlst = [];

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
    gF0 = @(x,theta0,theta1,theta2) (F0(x,theta0+dt,theta1,theta2,noise0) - F0(x,theta0-dt,theta1,theta2,noise0))/(2*dt);
    gF1 = @(x,theta0,theta1,theta2) (F1(x,theta0,theta1+dt,theta2,noise1) - F1(x,theta0,theta1-dt,theta2,noise1))/(2*dt);
    gF2 = @(x,theta0,theta1,theta2) (F2(x,theta0,theta1,theta2+dt,noise2) - F2(x,theta0,theta1,theta2-dt,noise2))/(2*dt);        

    % initialization
    Itheta0 = start(j,1); 
    Itheta1 = start(j,2); 
    Itheta2 = start(j,3);

    theta0 = Itheta0; 
    theta1 = Itheta1; 
    theta2 = Itheta2; 

    k = 1; % number of iteration

    interim_F0 = []; interim_F1 = []; interim_F2 = [];
    f_eval_F0 = 0; f_eval_F1 = 0; f_eval_F2 = 0;  
    total_iter = 0;
    
%     p_old = 0;
    
    while 1       
        p = randi([1 3],1);
        
%         while p == p_old
%             p = randi([1 3],1);
%         end
        
        % F0 block 
        if p == 1 
    
        m = 1;
        while m < 2
            [tau0,f_eval_F0] = back0(theta0,theta1,theta2,x0,x1,x2,noise0,noise1,noise2,case_no,f_eval_F0); % backtracking linesearch
%             tau0 = 0.01/(k+1)^2;
            
            theta0 = theta0 - tau0*gF0(x0,theta0,theta1,theta2); % parameter update          
            F0a = F0(x0,theta0,theta1,theta2,noise0);
            new_F0 = F0a; new_theta0 = theta0; f_eval_F0 = f_eval_F0 + 2;
            if m ~= 1
%                 delta0 = abs(new_F0-old_F0)/abs(old_F0);
                delta0 = norm(gF0(x0,theta0,theta1,theta2),2); f_eval_F0 = f_eval_F0 + 2;
                if (delta0<epsilon) || (m>max_iter_in) 
                    break;
                end
            end
            old_F0 = new_F0; old_theta0 = new_theta0; 
            m = m + 1;        
        end

        F0a = F0(x0,theta0,theta1,theta2,noise0);
        F1a = F1(x1,theta0,theta1,theta2,noise1);
        F2a = F2(x2,theta0,theta1,theta2,noise2);

        interim_F0 = [interim_F0 F0a];
        interim_F1 = [interim_F1 F1a];        
        interim_F2 = [interim_F2 F2a];

        % F1 block 
        elseif p == 2
                     
        m = 1;
        while m < 2
            [tau1,f_eval_F1] = back1(theta0,theta1,theta2,x0,x1,x2,noise0,noise1,noise2,case_no,f_eval_F1);
%             tau1 = 0.01/(k+1)^2;
            
            theta1 = theta1 - tau1*gF1(x1,theta0,theta1,theta2);           
            F1a = F1(x1,theta0,theta1,theta2,noise1); 
            new_F1 = F1a; new_theta1 = theta1; f_eval_F1 = f_eval_F1 + 2; 
            if m ~= 1
%                 delta1 = abs(new_F1-old_F1)/abs(old_F1);
                delta1 = norm(gF1(x1,theta0,theta1,theta2),2); f_eval_F1 = f_eval_F1 + 2;
                if (delta1<epsilon) || (m>max_iter_in) 
                    break;
                end
            end
            old_F1 = new_F1; old_theta1 = new_theta1; 
            m = m + 1;  
            total_iter = total_iter + 1;
        end

        F0a = F0(x0,theta0,theta1,theta2,noise0);
        F1a = F1(x1,theta0,theta1,theta2,noise1);
        F2a = F2(x2,theta0,theta1,theta2,noise2);

        interim_F0 = [interim_F0 F0a];
        interim_F1 = [interim_F1 F1a];        
        interim_F2 = [interim_F2 F2a];

        % F2 block
        elseif p == 3

        m = 1;
        while m < 2
            [tau2,f_eval_F2] = back2(theta0,theta1,theta2,x0,x1,x2,noise0,noise1,noise2,case_no,f_eval_F2);
%             tau2 = 0.01/(k+1)^2;
            
            theta2 = theta2 - tau2*gF2(x2,theta0,theta1,theta2);
            F2a = F2(x2,theta0,theta1,theta2,noise2);
            new_F2 = F2a; new_theta2 = theta2; f_eval_F2 = f_eval_F2 + 2;
            if m ~= 1
%                 delta2 = abs(new_F2-old_F2)/abs(old_F2);
                delta2 = norm(gF2(x2,theta0,theta1,theta2),2); f_eval_F2 = f_eval_F2 + 2;
                if (delta2<epsilon) || (m>max_iter_in) 
                    break;
                end
            end
            old_F2 = new_F2; old_theta2 = new_theta2; 
            m = m + 1;
            total_iter = total_iter + 1;
        end

        F0a = F0(x0,theta0,theta1,theta2,noise0);
        F1a = F1(x1,theta0,theta1,theta2,noise1);
        F2a = F2(x2,theta0,theta1,theta2,noise2);

        interim_F0 = [interim_F0 F0a];
        interim_F1 = [interim_F1 F1a];        
        interim_F2 = [interim_F2 F2a];

        end 

        % training loss
        Fta = Ft(x_init,theta0,theta1,theta2,noise_init);
        F0a = F0(x0,theta0,theta1,theta2,noise0); f_eval_F0 = f_eval_F0 + 1;
        F1a = F1(x1,theta0,theta1,theta2,noise1); f_eval_F1 = f_eval_F1 + 1;
        F2a = F2(x2,theta0,theta1,theta2,noise2); f_eval_F2 = f_eval_F2 + 1;

        % testing loss
        Ftb = Ft(x_test_init,theta0,theta1,theta2,noise_test_init);
        F0b = F0(x0_test,theta0,theta1,theta2,noise0_test); f_eval_F0 = f_eval_F0 + 1;
        F1b = F1(x1_test,theta0,theta1,theta2,noise1_test); f_eval_F1 = f_eval_F1 + 1;
        F2b = F2(x2_test,theta0,theta1,theta2,noise2_test); f_eval_F2 = f_eval_F2 + 1;

        % convergence criteria    
        new_F0c = F0a; new_F1c = F1a; new_F2c = F2a;
        new_theta0 = theta0; new_theta1 = theta1; new_theta2 = theta2; 

        if k ~= 1
%             eta0 = abs(new_F0c-old_F0c)/abs(old_F0c);
%             eta1 = abs(new_F1c-old_F1c)/abs(old_F1c);
%             eta2 = abs(new_F2c-old_F2c)/abs(old_F2c);

            eta0 = norm(gF0(x0,theta0,theta1,theta2),2); f_eval_F0 = f_eval_F0 + 2;
            eta1 = norm(gF1(x1,theta0,theta1,theta2),2); f_eval_F1 = f_eval_F1 + 2;
            eta2 = norm(gF2(x2,theta0,theta1,theta2),2); f_eval_F2 = f_eval_F2 + 2;     
            
            eta = max(eta0,max(eta1,eta2));

            if (eta<epsilon) || (k>max_iter_out) 
                break;
            end
        end
        old_F0c = new_F0c; old_F1c = new_F1c; old_F2c = new_F2c;
        old_theta0 = new_theta0; old_theta1 = new_theta1; old_theta2 = new_theta2; 
        k = k + 1;
%         p_old = p;
        total_iter = total_iter + 1; 
    end
    tempT = [theta0 theta1 theta2];
    tempF = [F0a F1a F2a Fta];  
    tempF_test = [F0b F1b F2b Ftb];  
    tempF_eval = [f_eval_F0 f_eval_F1 f_eval_F2];
    start_theta_rslt = [start_theta_rslt; tempT];
    start_F_rslt = [start_F_rslt; tempF];
    start_F_test_rslt = [start_F_test_rslt; tempF_test];
    start_F_eval_rlst = [start_F_eval_rlst; tempF_eval];
    start_total_iter_rlst = [start_total_iter_rlst; total_iter];
    start_outer_iter_rlst = [start_outer_iter_rlst; k];  
    
    disp(['= result: ' num2str(theta0) ', ' num2str(theta1) ', ' num2str(theta2) ]) 

    F0_array{iTry,1} = interim_F0;
    F1_array{iTry,1} = interim_F1;
    F2_array{iTry,1} = interim_F2;

    end
    theta_mean_store = [theta_mean_store; mean(start_theta_rslt,1)];
    theta_std_store = [theta_std_store; std(start_theta_rslt,1)];
    F_mean_store = [F_mean_store; mean(start_F_rslt,1)];
    F_std_store = [F_std_store; std(start_F_rslt,1)];
    F_test_mean_store = [F_test_mean_store; mean(start_F_test_rslt,1)];
    F_test_std_store = [F_test_std_store; std(start_F_test_rslt,1)]; 
    F_eval_mean_store = [F_eval_mean_store; mean(start_F_eval_rlst,1)];
    F_eval_std_store = [F_eval_std_store; std(start_F_eval_rlst,1)];
    total_iter_mean_store = [total_iter_mean_store; mean(start_total_iter_rlst,1)];
    total_iter_std_store = [total_iter_std_store; std(start_total_iter_rlst,1)];
    outer_iter_mean_store = [outer_iter_mean_store; mean(start_outer_iter_rlst,1)];
    outer_iter_std_store = [outer_iter_mean_store; std(start_outer_iter_rlst,1)];
    
end       

theta_mean_store
theta_std_store
F_mean_store
F_std_store
F_test_mean_store
F_test_std_store

F_eval_mean_store
F_eval_std_store % function evaluations 
total_iter_mean_store
total_iter_std_store % outer + inner iterations
outer_iter_mean_store
outer_iter_std_store % outer interations
time = toc

% select the result whose traning loss is minimum 
% [~,idx] = min(F_mean_store,[],1);
% idx = idx(1);
% result = [theta_mean_store(idx,:) theta_std_store(idx,:) ...
%             F_test_mean_store(idx,:) F_test_std_store(idx,:) ...
%             F_eval_mean_store(idx,:) F_eval_std_store(idx,:) ...
%             total_iter_mean_store(idx,:) total_iter_std_store(idx,:) ...
%             outer_iter_mean_store(idx,:) outer_iter_std_store(idx,:) time 0 0 0];

% #########################################################################
% 3. Save file
% #########################################################################

save(['/results/EX' num2str(case_no) '_M_BC_VS.mat'])
    
    