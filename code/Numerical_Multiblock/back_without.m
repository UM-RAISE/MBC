function [tau,f_eval_F] = back_without(theta,x0,noise0,case_no,f_eval_F)

dt = 10^(-8); 

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

theta0 = theta(1);
theta1 = theta(2);
theta2 = theta(3);

% problem 
if case_no == 1
    fpt = @(x,theta0,theta1,theta2,noiset) Ia(x).*exp(x./10).*sin(x) + Ib(x).*exp(x./10).*cos(x) + Ic(x).*exp(x./10).*cos(x) + Id(x).*exp(x./10).*cos(x) + Ie(x).*exp(x./10).*sin(x) + noiset;
    fp0 = @(x,theta0,theta1,theta2,noise0) Ia(x).*exp(x./10).*sin(x) + Ib(x).*exp(x./10).*cos(x) + noise0; 
    fp1 = @(x,theta0,theta1,theta2,noise1)                             Ib(x).*exp(x./10).*cos(x) + Ic(x).*exp(x./10).*cos(x) + Id(x).*exp(x./10).*cos(x) + noise1;
    fp2 = @(x,theta0,theta1,theta2,noise2)                                                                                     Id(x).*exp(x./10).*cos(x) + Ie(x).*exp(x./10).*sin(x) + noise2; 

    fct = @(x,theta0,theta1,theta2) exp(x./10).*sin(x) + exp(x./10).*cos(x) + exp(x./10).*sin(x)   - (1/2)*(theta0+1).*sqrt(x.^2-x+1) - (1/10).*(theta0+1)*(theta1-10).*sqrt(x.^2-x+1) - (1/10).*(theta1-10).*sqrt(x.^2-x+1) - (1/10).*(theta1-10)*(theta2-5).*sqrt(x.^2-x+1) - (1/20).*(theta2-5).*sqrt(x.^2-x+1);
    fc0 = @(x,theta0,theta1,theta2) exp(x./10).*sin(x) + exp(x./10).*cos(x)                        - (1/2)*(theta0+1).*sqrt(x.^2-x+1) - (1/10).*(theta0+1)*(theta1-10).*sqrt(x.^2-x+1); 
    fc1 = @(x,theta0,theta1,theta2)                      exp(x./10).*cos(x)                                                           - (1/10).*(theta0+1)*(theta1-10).*sqrt(x.^2-x+1) - (1/10).*(theta1-10).*sqrt(x.^2-x+1) - (1/10).*(theta1-10)*(theta2-5).*sqrt(x.^2-x+1);
    fc2 = @(x,theta0,theta1,theta2)                      exp(x./10).*cos(x) + exp(x./10).*sin(x)                                                                                                                             - (1/10).*(theta1-10)*(theta2-5).*sqrt(x.^2-x+1) - (1/20).*(theta2-5).*sqrt(x.^2-x+1);            

elseif case_no == 2
    fpt = @(x,theta0,theta1,theta2,noiset) Ia(x).*exp(x./10).*sin(x) + Ib(x).*exp(x./10).*cos(x) + Ic(x).*exp(x./10).*cos(x) + Id(x).*exp(x./10).*cos(x) + Ie(x).*exp(x./10).*sin(x) + noiset;
    fp0 = @(x,theta0,theta1,theta2,noise0) Ia(x).*exp(x./10).*sin(x) + Ib(x).*exp(x./10).*cos(x) + noise0; 
    fp1 = @(x,theta0,theta1,theta2,noise1)                             Ib(x).*exp(x./10).*cos(x) + Ic(x).*exp(x./10).*cos(x) + Id(x).*exp(x./10).*cos(x) + noise1;
    fp2 = @(x,theta0,theta1,theta2,noise2)                                                                                     Id(x).*exp(x./10).*cos(x) + Ie(x).*exp(x./10).*sin(x) + noise2; 

    fct = @(x,theta0,theta1,theta2) exp(x./10).*sin(x) + exp(x./10).*cos(x) + exp(x./10).*sin(x)   - (2)*(theta0+1).*(cos(theta0.*x./10)) - (1/2).*(theta0+1)*(theta1-10).*(cos(theta1.*x./20)) - (1/2)*(theta1-10).*(cos(theta1.*x./10)) - (1/2).*(theta1-10)*(theta2-5).*(cos(theta1.*x./20)) - (1/2).*(theta2-5).*(cos(theta2.*x./10));
    fc0 = @(x,theta0,theta1,theta2) exp(x./10).*sin(x) + exp(x./10).*cos(x)                        - (2)*(theta0+1).*(cos(theta0.*x./10)) - (1/2).*(theta0+1)*(theta1-10).*(cos(theta1.*x./20)); 
    fc1 = @(x,theta0,theta1,theta2)                      exp(x./10).*cos(x)                                                               - (1/2).*(theta0+1)*(theta1-10).*(cos(theta1.*x./20)) - (1/2)*(theta1-10).*(cos(theta1.*x./10)) - (1/2).*(theta1-10)*(theta2-5).*(cos(theta1.*x./20));
    fc2 = @(x,theta0,theta1,theta2)                      exp(x./10).*cos(x) + exp(x./10).*sin(x)                                                                                                                                          - (1/2).*(theta1-10)*(theta2-5).*(cos(theta1.*x./20)) - (1/2).*(theta2-5).*(cos(theta2.*x./10));

elseif case_no == 3
    fpt = @(x,theta0,theta1,theta2,noiset) exp(x./10) + I(x).*exp(x./10).*sin(x) + (1-I(x)).*exp(x./10).*cos(x) + noiset;
    fp0 = @(x,theta0,theta1,theta2,noise0) exp(x./10) + I(x).*exp(x./10).*sin(x) + (1-I(x)).*exp(x./10).*cos(x) + noise0; 
    fp1 = @(x,theta0,theta1,theta2,noise1) exp(x./10) + exp(x./10).*sin(x) + noise1;
    fp2 = @(x,theta0,theta1,theta2,noise2) exp(x./10) + exp(x./10).*cos(x) + noise2; 

    fct = @(x,theta0,theta1,theta2) exp(x./10) + exp(x./10).*sin(x) + exp(x./10).*cos(x) - 5.*sqrt(theta0.^2-theta0+1).*(sin(theta0.*x./10)+cos(theta0.*x./10)) - (1/2).*sqrt(theta1.^2-theta1+1).*(sin(theta1.*x./5)+cos(theta1.*x./5)) - (1/2).*sqrt(theta2.^2-theta2+1).*(sin(theta2.*x./10)+cos(theta2.*x./10));
    fc0 = @(x,theta0,theta1,theta2) exp(x./10) + exp(x./10).*sin(x) + exp(x./10).*cos(x) - 5.*sqrt(theta0.^2-theta0+1).*(sin(theta0.*x./10)+cos(theta0.*x./10)) - (1/2).*sqrt(theta1.^2-theta1+1).*(sin(theta1.*x./5)+cos(theta1.*x./5)) - (1/2).*sqrt(theta2.^2-theta2+1).*(sin(theta2.*x./10)+cos(theta2.*x./10));
    fc1 = @(x,theta0,theta1,theta2) exp(x./10) + exp(x./10).*sin(x) - 5.*sqrt(theta0.^2-theta0+1).*(sin(theta0.*x./10)+cos(theta0.*x./10)) - (1/2).*sqrt(theta1.^2-theta1+1).*(sin(theta1.*x./5)+cos(theta1.*x./5)); 
    fc2 = @(x,theta0,theta1,theta2) exp(x./10) + exp(x./10).*cos(x) - 5.*sqrt(theta0.^2-theta0+1).*(sin(theta0.*x./10)+cos(theta0.*x./10)) - (1/2).*sqrt(theta2.^2-theta2+1).*(sin(theta2.*x./10)+cos(theta2.*x./10));                
end   

% objective function
Ft = @(x,theta0,theta1,theta2,noiset) (1/length(x))*sum((fct(x,theta0,theta1,theta2) - fpt(x,Ttheta0,Ttheta1,Ttheta2,noiset)).^2);
F0 = @(x,theta0,theta1,theta2,noise0) (1/length(x))*sum((fc0(x,theta0,theta1,theta2) - fp0(x,Ttheta0,Ttheta1,Ttheta2,noise0)).^2);
F1 = @(x,theta0,theta1,theta2,noise1) (1/length(x))*sum((fc1(x,theta0,theta1,theta2) - fp1(x,Ttheta0,Ttheta1,Ttheta2,noise1)).^2);
F2 = @(x,theta0,theta1,theta2,noise2) (1/length(x))*sum((fc2(x,theta0,theta1,theta2) - fp2(x,Ttheta0,Ttheta1,Ttheta2,noise2)).^2);

% gradient
% gF0 = @(x,theta0,theta1,theta2) (F0(x,theta0+dt,theta1,theta2,noise0) - F0(x,theta0-dt,theta1,theta2,noise0))/(2*dt);
% gF1 = @(x,theta0,theta1) (F1(x,theta0,theta1+dt,noise1) - F1(x,theta0,theta1-dt,noise1))/(2*dt);
% gF2 = @(x,theta0,theta2) (F2(x,theta0,theta2+dt,noise2) - F2(x,theta0,theta2-dt,noise2))/(2*dt);      

d1 = @(x,theta0,theta1,theta2,noise0) (Ft(x,theta0+dt,theta1,theta2,noise0) - Ft(x,theta0-dt,theta1,theta2,noise0))/(2*dt); 
d2 = @(x,theta0,theta1,theta2,noise0) (Ft(x,theta0,theta1+dt,theta2,noise0) - Ft(x,theta0,theta1-dt,theta2,noise0))/(2*dt);
d3 = @(x,theta0,theta1,theta2,noise0) (Ft(x,theta0,theta1,theta2+dt,noise0) - Ft(x,theta0,theta1,theta2-dt,noise0))/(2*dt);

gF1 = [d1(x0,theta0,theta1,theta2,noise0);d2(x0,theta0,theta1,theta2,noise0);d3(x0,theta0,theta1,theta2,noise0)];
f_eval_F = f_eval_F + 7; % one function eval + six function evals with h

% backtracking linesearch 
alpha_bar = 1;
rho = 0.5;
c1 = 10^(-4);
alpha = alpha_bar;
while (Ft(x0,theta0-alpha*gF1(1),theta1-alpha*gF1(2),theta2-alpha*gF1(3),noise0) > Ft(x0,theta0,theta1,theta2,noise0) - c1*alpha*norm(gF1,2)^2)
    alpha = alpha*rho;
    f_eval_F = f_eval_F + 1;
end    
tau = alpha;

end