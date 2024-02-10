clc
clear
close all
format long

syms x
%% Initializing constants
a = 1;  % goal is to approximate f'(a)
numPts = 4;
numH = 21;

%% Constructing h and x values

% constructing uniform increments h
h = ones(1, numH);
h = 10*h;
for i = 1:numH
    h(i) = h(i)^(1-i);
end

% constructing x values using different increments h
X = zeros(numPts, numH);
X(1,:) = a*ones(1, numH);  % a = entries of 1st row of X
for j = 1:numH
    for i = 2:numPts
        X(i,j) = X(1,j) - (i-1)*h(j);
    end
end

%% f, estimates for f'(a) and relative errors
f = exp(X);  % functions values
exact = exp(1);  % exact e
D = zeros(1, numH);

% estimates D(h)
for j = 1:numH
    num = f(4,j) - f(3,j) - 5*f(2,j) + 5*f(1,j);
    den = 4*h(j);
    D(j) = num/den;  % approximating formula
end

disp('estimates D(h) = ')
disp(vpa(D, 6)')

% relative errors for D(h)
relError = zeros(1, numH);
for j = 1:numH
    numer = abs(exact - D(j));
    denom = abs(exact);
    relError(j) = numer/denom;
end

disp('rel Error (h) = ')
disp(vpa(relError*100, 6)')

%% solving new estimates using stepsize h/2
newh = h/2;

% constructing x values using h/2
newX = zeros(numPts, numH);
newX(1,:) = a*ones(1, numH);  % a = entries of 1st row of X
for j = 1:numH
    for i = 2:numPts
        newX(i,j) = newX(1,j) - (i-1)*newh(j);
    end
end

newf = exp(newX);
D_2 = zeros(1, numH);
% estimates D(h/2)
for j = 1:numH
    num = newf(4,j) - newf(3,j) - 5*newf(2,j) + 5*newf(1,j);
    den = 4*newh(j);
    D_2(j) = num/den;  % approximating formula
end


%% New Estimates using Richardson Extrapolation
% formula obtained:
% 4D(h/2) - D(h))/3
RE_est = zeros(1, numH);
for j = 1:numH
    RE_est(j) = (4*D_2(j) - D(j))/3;
end

disp('estimates Re_est = ')
disp(vpa(RE_est, 6)')

% relative errors for RE_est
newrelError = zeros(1, numH);
for j = 1:numH
    numer = abs(exact - RE_est(j));
    denom = abs(exact);
    newrelError(j) = numer/denom;
end

disp('rel Error (RE_est) = ')
disp(vpa(newrelError*100, 6)')