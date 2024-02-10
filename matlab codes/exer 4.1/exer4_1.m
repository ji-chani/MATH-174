% 4-point backward formula of order three
clc
clear
format long

syms x
%%  Initializing Constants and Functions
f = exp(x);  % function to approximate
a = 1;  % where f'(a) is the value to approximate
numh = 21;  % no. of increments
numPts = 4;  % no. of data points

%% Icrements h and x values

% constructing uniform increments h
h = ones(1, numh);
h = 10*h;
for i = 1:numh
    h(i) = h(i)^(1-i);
end

% constructing x's using different h
% the output is a numPts x numh matrix
% each column is a group of x's for a specific h
% the h value for each column in X is in 
% the same indexed column of the h array
X = zeros(numPts, numh);
X(1,:) = a*ones(1, numh);  % entries in first row of X = a
for j = 1:numh
    for i = 2:numPts  % starts iteration on 2nd row
        X(i, j) = X(1,j) - (i-1)*(h(j));  % (x-h), (x-2h), (x-3h)
    end
end

F = exp(X);  % f(x) = e^x

%%  Calculation of f'(a) and relative error

% f'(a)
Fprime = zeros(1, numh);
for j = 1:numh
    num = 11*F(1,j) - 18*F(2,j) + 9*F(3,j) - 2*F(4,j);
    den =6*h(j);
    Fprime(j) = num/den;  % approximating formula
end
disp("f'(1) = ")
disp(Fprime)

% relative error
relError = zeros(1, numh);
exact = double(subs(f, 1));
for j = 1:numh
    numer = exact - Fprime(j);
    relError(j) = abs(numer/exact);
end

disp('relative errors (in %) =')
disp(vpa(relError*100, 6)')
