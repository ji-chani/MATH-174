clc
clear
close all

syms x real
syms S F
%% Abscissas
X = 0:0.1:1;
numPts = length(X);  % n+1 (no. of data points)

% calculation of Y
mn = 0.5204;
sd = 0.2196;

for i = 1:numPts
    denom = sd * sqrt(2*pi);
    num = (X(i) - mn)^2;
    den = 2 * sd^2;
    n = num/den;
    F(i) = (1/denom) * exp(-n);
end

%% Computation of constants a and increment h
a = F';  % constant a
h = zeros(1, numPts-1);  % increments in each interval of abscissas
for i = 1:length(h)
    h(i) = X(i+1)-X(i);
end

%% Computation of constants c using matrices
A = zeros(numPts, numPts);  % square matrix A
B = zeros(numPts, 1);  % column vector b

for k = 2:numPts-1
    % supplying middle rows of A
    A(k, k-1) = h(k-1);
    A(k, k) = 2*(h(k-1)+h(k));
    A(k, k+1) = h(k);
       
    % supplying middle rows of B
    left = (3*(a(k+1)-a(k)))/h(k);
    right = (3*(a(k)-a(k-1)))/h(k-1);
    B(k) = left - right;
end

% imposing free conditions (not-a-knot)
tb_rows = zeros(numPts, numPts);  % top and bottom rows of A
for k = linspace(2, numPts-1, 2)  % iterates over 2 and n-1 only
    tb_rows(k, k-1) = h(k);
    tb_rows(k, k) = -(h(k-1)+h(k));
    tb_rows(k, k+1) = h(k-1);
end
A(1,:) = tb_rows(2,:);  % 1st row of A becomes 2nd row of tb_rows
A(numPts,:) = tb_rows(numPts-1,:);  % n+1th row of A becomes nth row of tb_rows

% computes for c (inv(A)*B)
c = A\B;

%% Computation of constants b and d using c

d = zeros(numPts-1, 1);  % column vector of size n-1
b = zeros(numPts-1, 1);  % column vector of size n-1
for k = 1:numPts-1
    % constant d
    d(k) = (c(k+1)-c(k))/(3*h(k));
    
    % constant b
    l = (a(k+1)-a(k))/h(k);
    r = (2*c(k)+c(k+1))/3 * h(k);
    b(k) = l-r;
end

%% Constructing Piecewise Cubic Spline
for i = 1:numPts-1
    S(i) = a(i) + b(i)*(x-X(i)) + c(i)*(x-X(i))^2 + d(i)*(x-X(i))^3;
end
S = S';
disp(vpa(expand(S), 4))

%% Plotting cubic splines
s = piecewise( ...
    (X(1)<=x)&(x<=X(2)), S(1), ...
    (X(2)<=x)&(x<=X(3)), S(2), ...
    (X(3)<=x)&(x<=X(4)), S(3), ...
    (X(4)<=x)&(x<=X(5)), S(4), ...
    (X(5)<=x)&(x<=X(6)), S(5), ...
    (X(6)<=x)&(x<=X(7)), S(6), ...
    (X(7)<=x)&(x<=X(8)), S(7), ...
    (X(8)<=x)&(x<=X(9)), S(8), ...
    (X(9)<=x)&(x<=X(10)), S(9), ...
    (X(10)<=x)&(x<=X(11)), S(10));

figure
hold on
fplot(s, 'r')
scatter(X, F, 'filled', 'bo')