clc
clear
format long

syms x
%% Initialization of constants

% integration bounds
a = 0;
b = 4;

% no. of data points
N1 = 100;
N2 = 25000;

% no. of trials
trials = 5;
t = 1;  % for 1st iteration

% function and estimates
f = x*sin(x^3);
est = zeros(2, trials);  % placeholder for estimates

%% Generating random numbers and Estimates
while t <= trials
    % uniformly distributed random numbers over (0,1)
    x1 = rand(N1, 1);  % random number-filled N1x1 matrix
    x2 = rand(N2, 1);

    % translating bounds
    x1 = (b-a).*x1 + a;  % r = (b-a)x + a
    x2 = (b-a).*x2 + a;
    
    % function values
    F1 = x1.*sin(x1.^3);
    F2 = x2.*sin(x2.^3);
    
    % estimates
    est(1, t) = (b-a) * mean(F1);
    est(2, t) = (b-a) * mean(F2);
    
    % counter
    t = t+1;
end
disp(est)