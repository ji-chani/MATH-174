clc
clear
close all

syms x real
syms S
%% Function and Necessary Constants
f = 1/(sqrt(2*pi))*exp(-(x^2/2));
interval =[-3, 3];

%% Finding increment h
ub = 0.0156;  % upper bound L^infty norm
f4Prime = diff(f, 4);
c = 5/384;

% plot of f4Prime
figure
fplot(abs(f4Prime), interval, 'r')
xlabel("$x$", 'Interpreter', 'Latex')
ylabel("$|f^{(4)}(x)|$", 'Interpreter', 'Latex')

% maximum value of abs(f4Prime) occurs at x=0
maxf4Prime = double(abs(subs(f4Prime, 0)));

% calculation of h (uniform increment)
h = abs(nthroot((ub)/(c*maxf4Prime), 4));
disp('h = ')
disp(h)

h = 1.00;  % after rounding down to 2 decimal places
disp('h = ')
disp(h)

%% Abscissas
X = min(interval):h:max(interval);  % construct an array with stepsize h and endpoints min(), max()
Y = 1/(sqrt(2*pi))*exp(-(X.^2/2));
numPts = length(X);  % n+1 (no. of abscissas)
disp('x =')
disp(X)
disp('y = ')
disp(Y')

% computation for constant a
a = Y;
%% Solving for constant c
A = zeros(numPts, numPts);  % matrix A of eqn (12)
B = zeros(numPts, 1);  % matrix b of eqn (12)
for k = 2:numPts-1  % iterates over 2 to n
    % supplying middle rows of A
    A(k, k-1) = h;
    A(k, k) = 4*h;
    A(k, k+1) = h;

    % supplying middle rows of b
    B(k) = (3/h) * (a(k+1) - 2*a(k) + a(k-1));
end

% imposing free conditions (clamped)
fPrime = diff(f, 1);  % first derivative of f
fPrimea = double(subs(fPrime, -3));  %fPrime(-3)
fPrimeb = double(subs(fPrime, 3));  %fPrime(3)
A(1, 1) = 2*h;
A(1, 2) = h;
A(numPts, numPts-1) = h;
A(numPts, numPts) = 2*h;
B(1) = (3/h)*(a(2)-a(1)) - (3*fPrimea);
B(numPts) = (3*fPrimeb) - (3/h) * (a(numPts)-a(numPts-1));

% computes for c
c = A\B;
%% Solving for constants b and d
b = zeros(numPts-1, 1);
d = zeros(numPts-1, 1);
for k = 1:numPts-1
    % constant d
    d(k) = (c(k+1)-c(k))/(3*h);

    % constant b
    l = (a(k+1)-a(k))/h;
    r = (2*c(k)+c(k+1))/3 * h;
    b(k) = l-r;
end
disp('b =')
disp(b)
disp('d =')
disp(d)

%% Constructing piecewise cubic splines
for i = 1:numPts-1
    S(i) = a(i) + b(i)*(x-X(i)) + c(i)*(x-X(i))^2 + d(i)*(x-X(i))^3;
end

disp(vpa(expand(S'), 4))

%% Plotting S with f
s = piecewise( ...
    (X(1)<=x)&(x<=X(2)), S(1), ...
    (X(2)<=x)&(x<=X(3)), S(2), ...
    (X(3)<=x)&(x<=X(4)), S(3), ...
    (X(4)<=x)&(x<=X(5)), S(4), ...
    (X(5)<=x)&(x<=X(6)), S(5), ...
    (X(6)<=x)&(x<=X(7)), S(6));

figure
hold on
fplot(s, 'r')
fplot(f, interval, 'k')
% scatter(X, Y, 'filled', 'k')  % scatter plot of data points (optional)
grid
legend('clamped cubic spline $s(x)$', '$f(x)$', ...
    'Interpreter', 'Latex', 'Fontsize', 14)
xlabel('$x$', ...
    'Interpreter', 'Latex', 'Fontsize', 14)
ylabel('$y$', ...
    'Interpreter', 'Latex', 'Fontsize', 14)
xlim([-3.1, 3.1])