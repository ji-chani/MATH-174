% Lagrange Form of the AIP
clear
close all

syms x % Declare x as a symbolic expression
%% Initialize Abscissas and Function Values
T = [0 10 25 50 100];
D = [0 8 4251 20677 357];

% T = [1 2 3 4 5 6];
% D = [1 2 3 8 13];

numPts = size(T, 2);  % represents n+1 (no. of data points)
interpInt = [min(T), max(T)];  % interpolating interval

%% Compute the AIP
P = 0;  % since P is a sum of terms, we initialize the sum to be 0
for k = 1:numPts
    lk = 1;  % Lagrange basis polynomial (product so we initialize value as 1)

    % Compute for Lagrange basis polynomial
    for i = 1:numPts
        if i ~= k
            numer = x - T(i);
            denom = T(k) - T(i);
            lk = lk * (numer/denom);
        end
    end
    P = P + D(k) * lk;
end

%% Postprocessing

% display symbolic or vpa representation
P = simplify(P);
disp('The AIP relative to the data points is')
pretty(P);  % writes the polynomial pretty-wise
disp('or')
disp(vpa(expand(P), 4))  % vpa returns values in decimal form (4 sig figs)
% expand(P) distributes the denominators of P to every ter

% display graph
figure
fplot(P, interpInt)  % if plot is symbolic
hold on
scatter(T, D, 'filled')  % 'filled' fills marker with color
title('Plot of the data points and/or the function, and its interpolant')
lgnd = legend('Interpolant $P$', 'Data Points', 'Interpreter', 'Latex');
xlabel('t')
ylabel('number of infected individuals')
set(lgnd, 'FontSize', 14)

% create a table of values
