clc
clear
close all

syms x real
syms S
%% Abscissas
X = [72.00 81.67 90.33 90.83 91.67 92.50 93.33 95.00 97.50 97.67];
Y = [80.57 85.20 83.57 87.37 95.68 99.71 95.79 98.00 94.77 94.40
];
numPts = length(X);  % n+1

%% Computation of increment h and constant a
h = zeros(1, numPts-1);  % h(i) = X(i+1)-X(i) for i = 1 to n
for i = 1:length(h)
    h(i) = X(i+1)-X(i);
end
a = Y';

%% Computation of constant c using matrices
A = zeros(numPts, numPts);  % square matrix A (n+1 by n+1)
B = zeros(numPts, 1);  % column vector B (n+1 by 1)

for k = 2:numPts-1  % iterates over 2 to n
    % supplying middle rows of A
    A(k, k-1) = h(k-1);
    A(k, k) = 2*(h(k-1)+h(k));
    A(k, k+1) = h(k);

    % suppying middle rows of B
    left = (3*(a(k+1)-a(k)))/h(k);
    right = (3*(a(k)-a(k-1)))/h(k-1);
    B(k) = left-right;
end

% imposing free conditions (natural)
% s''(X(1)) = s''(X(n+1)) = 0
% => 2c(i) = 0 and 2(c(n+1)) = 0
for k = [1 numPts]
    A(k, k) = 2;
end

% computes for c (c = inv(A)*B)
c = A\B;

%% Computes for constants b and d
b = zeros(numPts-1, 1);  % column vector of size n
d = zeros(numPts-1, 1);  % column vector of size n
for k = 1:numPts-1
    % constant d
    d(k) = (c(k+1)-c(k))/(3*h(k));

    % constant b
    l = (a(k+1)-a(k))/h(k);
    r = (2*c(k)+c(k+1))/3 * h(k);
    b(k) = l-r;
end

%% Constructing Piecewise Cublic Splines
for i = 1:numPts-1
    S(i) = a(i) + b(i)*(x-X(i)) + c(i)*(x-X(i))^2 + d(i)*(x-X(i))^3;
end

disp(vpa(expand(S'), 4))

%% Plotting Cubic Splines
s = piecewise( ...
    (X(1)<=x)&(x<=X(2)), S(1), ...
    (X(2)<=x)&(x<=X(3)), S(2), ...
    (X(3)<=x)&(x<=X(4)), S(3), ...
    (X(4)<=x)&(x<=X(5)), S(4), ...
    (X(5)<=x)&(x<=X(6)), S(5), ...
    (X(6)<=x)&(x<=X(7)), S(6), ...
    (X(7)<=x)&(x<=X(8)), S(7), ...
    (X(8)<=x)&(x<=X(9)), S(8), ...
    (X(9)<=x)&(x<=X(10)), S(9));

figure
hold on
fplot(s, 'r')
scatter(X, Y, 'filled', 'ko')
grid
legend('natural cubic spline $s(x)$', 'data points', ...
    'Interpreter', 'Latex', 'Fontsize', 14)
xlabel('$x$', ...
    'Interpreter', 'Latex', 'Fontsize', 14)
ylabel('$y$', ...
    'Interpreter', 'Latex', 'Fontsize', 14)
xlim([70, 100])
ylim([72, 102])

%% Substituting Estimate
est = 92.80;
pred_rating = double(subs(S(6), est));
disp(pred_rating)