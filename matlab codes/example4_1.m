clc
clear
close all

syms F S
syms x real

%% 
numNodes = 11;
T = linspace(0, 1, numNodes);
mn = 0.5204;
sd = 0.2196;

for i = 1:numNodes
    denom = sd * sqrt(2*pi);
    num = (T(i) - mn)^2;
    den = 2 * sd^2;
    n = num/den;
    F(i) = (1/denom) * exp(-n);
end

% linear polynomials S_i
for i = 1:numNodes-1
    a = F(i);
    b = (F(i+1)- F(i))/(T(i+1)-T(i));
    S(i) = a + b*(x - T(i));
end

disp(vpa(expand(S'), 4))

Est = 0;
for i =6:9
    new = double(int(S(i+1), T(i+1), T(i+2)));
    Est = Est + new;
end

disp(Est)

%%
% function func = P(z)  % returns func

%     ms = 0.5204;  % mean score
%     sd = 0.2196;  % standard deviation
%     denom = sd * sqrt(2 * pi);
%     n = (z - ms)^2 / (2*sd^2);
%     func = (1/denom) * exp(-n);
% end
