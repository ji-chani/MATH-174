clc
clear
close all

syms S
syms x real
%% Abscissas
numNodes = 7;
T = linspace(-3, 3, numNodes);
F = zeros(1,numNodes);  % F(x)
for i = 1:numNodes
    den = sqrt(2*pi);
    n = (T(i)^2/2);
    F(i) = (1/den)*exp(-n);
end

%% Calculation of linear interpolants S_i
for i = 1:numNodes-1
    a = F(i);
    b = (F(i+1)-F(i))/(T(i+1)-T(i));
    S(i) = a + b*(x-T(i));
end
S = S';
disp(vpa(expand(S), 4))