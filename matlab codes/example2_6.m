% This code is for Example 2.6 (MATH 174 - Module 2)

clear
close all


syms x  % Declares x as a symbolic variable
%% Defining Abscissas
T = linspace(0, 0.5, 9);
D = exp(T.^2);

numPts = length(T);
interval = [min(T), max(T)];

%% Compute for AIP
P = 0;
for k = 1:numPts
    lk = 1;
    for i = 1:numPts
        if i ~= k
            numer = x - T(i);
            denom = T(k) - T(i);
            lk = lk * (numer/denom);
        end
    end
    P = P + D(k) * lk;
end

%% Display Function and its Plot

% function
P = simplify(P);
disp('The AIP relative to the data points is')
disp(vpa(expand(P), 4))

% plot
figure
fplot(P, interval, 'k')
hold on
scatter(T, D, 'r', 'filled')
grid

xlim([-0.1, 0.6])
ylim([0.9 1.35])

%% Integrating P over the interval [0, 1/2]
disp('and the definite integral from x = 0 to x = 0.5')
disp('of e^(x^2) is equal to')
double(int(P, 0, 0.5))