clear
close all

x = sym('x');  % declares the symbolic variable x
%% Initializes abscissas and function values

% closing price of ICT stocks in PH
% obtained from https://ph.investing.com/equities/intl-container-historical-data
% sep 1 - 182.90
% aug 1 - 183.90
% july 1 - 191.00
% june 1 - 210.20
% may 2 - 215.60
% april 1 - 226.60
% march 1 - 218.00
% feb 2 - 200.20
% jan 3 - 195.00

T = [3 33 60 91 122 152 182 213 244];  % time variable T
S = [182.90 183.90 191.00 210.20 215.60 226.60 218.00 200.20 195.00];  % stock prices S

numPts = size(T, 2);
interpInt = [1 246];

%% Compute the AIP

% since P is a sum of terms, we initialize P as
P = 0;
for i = 1:numPts
    % since the Lagrange basis polynomial is a product of terms
    % we initialize lk as
    lk = 1;

    % Computation of lk
    for j = 1:numPts
        if i ~= j
            numer = x - T(j);
            denom = T(i) - T(j);
            lk = lk * (numer/denom);
        end
    end
    P = P + S(i) * lk;
end

%% Postprocessing

% Displaying the polynomial
P = simplify(P);
disp('The AIP relative to the data point is')
pretty(P);  % displays P decently

disp('or')
disp(vpa(expand(P), 4))  % dislays P having decimal form coeffients w 4 sf


% Displaying the plots
figure
fplot(P, interpInt, 'k')  % creates a 2-D plot of P (y-values) over interInt (x-interval)
hold on
scatter(T, S, 'filled')  % scatter plot of T vs S
grid

title('Plot of data points and its interpolant')
lgnd = legend('Interpolant $P$', 'Data Points', 'Interpreter', 'Latex');
x_l = xlabel('time variable ($t$)', 'Interpreter', 'Latex');
y_l = ylabel('closing stock price (in pesos)', 'Interpreter', 'Latex');

for k = [lgnd x_l y_l]
    set(k, 'Fontsize', 14)
end

