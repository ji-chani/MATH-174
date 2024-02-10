clc
clear
close all

syms x L 
%% 
a = 2;
N = 10^6;
f = 1/(log(x));  % f(x) = 1/ln(x), log in MATLAB is natural log

nP5 = 5;

%% Obtaining Legendre Roots
L0 = 1;
L(1) = x;
L(2) = 3/2*x^2 - 1/2;
for n = 3:nP5
    left = ((2*n-1)/n)*x*L(n-1);
    right = ((n-1)/n)*L(n-2);
    L(n) = left - right;  % L_n = (2n-1)/n xL_(n-1) - (n-1/n)L(n-2)
end
L4 = L(4);  % degree 4 Legendre polynomial 
L5 = L(5);  % degree 5 Legendre polynomial

T = zeros(2,nP5);  % abscissas
T(1,:) = [roots(sym2poly(L4))', 0];  % deg 4 roots
T(2,:) = roots(sym2poly(L5))';  % deg 5 roots

%% Obtaining Weights

w = zeros(2, nP5);
L4prime = diff(L4);
L5prime = diff(L5);
for i = 1:nP5
    d1 = double(subs(L4prime, T(1,i)));
    d2 = double(subs(L5prime, T(2,i)));
    
    den1 = (1-T(1,i)^2)*d1^2;
    den2 = (1-T(2,i)^2)*d2^2;

    w(1,i) = 2/den1;
    w(2,i) = 2/den2;
end

%% Gauss-Legendre Quadrature Rules
c = (N-a)/2;
d = (a+N)/2;
f_t = (log(c*T + d)).^(-1);  %f(ct+d)

sum1 = 0;
for i = 1:nP5-1
    sum1 = sum1 + w(1,i)*f_t(1,i);
end

sum2 = 0;
for i = 1:nP5
    sum2 = sum2 + w(2,i)*f_t(2,i);
end
Est = zeros(2,1);
Est(1) = c*sum1;
Est(2) = c*sum2;
disp(Est)

%% Calculation of Relative Error
exact = 78498;
relError = abs((Est-exact)/exact);

disp(relError * 100)