clc
clear
close all
syms x real
syms S
%% Abscissas
T = [3 4 5 7 8 9 10 11];
T = 100*T;
numPts = length(T);  % no of data points
W = [0.024 0.035 0.046 0.058 0.067 0.083 0.097 0.111];

%% Solving for Cubic Splines
a = W;
h = zeros(1, numPts-1);  % increments in each interval
for i = 1:length(h)
    h(i) = T(i+1) - T(i);
end

% defining necessary matrices
A = zeros(numPts, numPts);  % square matrix A of size numPts
b = zeros(numPts, 1);  % column vector b of size numPts

% suppying middle rows of A and b
for k = 2:numPts-1
    % matrix A
    A(k, k-1) = h(k-1);
    A(k, k) = 2*(h(k-1) + h(k));
    A(k, k+1) = h(k);
    
    % vector b
    left = (3*(a(k+1)-a(k)))/h(k);
    right = (3*(a(k)-a(k-1)))/h(k-1);
    b(k) = left - right;
end

% imposing free conditions (not-a-knot)
% supplying top and bottom rows of A
tb_rows = zeros(2, numPts);
for i = linspace(2, numPts-1, 2) % iterates over 2 and numPts-1 only
    if i == 2  
        tb_rows(1, i-1) = h(i);
        tb_rows(1, i) = -(h(i-1)+h(i));
        tb_rows(1, i+1) = h(i-1);
    else
        tb_rows(2, i-1) = h(i);
        tb_rows(2, i) = -(h(i-1)+h(i));
        tb_rows(2, i+1) = h(i-1);
    end
end
A(1,:) = tb_rows(1,:);
A(numPts,:) = tb_rows(2,:);

% computing for c_i
c = A\b;

% computing b_i and d_i
b = zeros(numPts-1, 1);  % matrix of size nx1
d = zeros(numPts-1, 1);  % matrix of size nx1
for i = 1:numPts-1
    l = (a(i+1)-a(i))/h(i);
    r = (2*c(i)+c(i+1))/3 * h(i);
    b(i) = l-r;
    d(i) = (c(i+1)-c(i))/(3*h(i));
end

% constructing polynomial
for i = 1:numPts-1
    S(i) = a(i) + b(i)*(x-T(i)) + c(i)*(x-T(i))^2 + d(i)*(x-T(i))^3;
end
S = S';
disp(vpa(expand(S), 4))