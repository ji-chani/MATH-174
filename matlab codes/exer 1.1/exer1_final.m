clearvars -except data  % to clear all existing variables except imported data (historical data)
close all

t = sym('t');  % declares t as a symbolic variable

%% Item 1: Constructing AIP (P)
% Initializes abscissas and function values

T = [3 33 60 91 122 152 182 213 244];  % time variable t
S = [195 200.20 218 226.6 215.6 210.2 191 183.9 182.9];  % stock prices S

numPts = size(T, 2);
interval = [1, 246];

% Compute for AIP

% since P is a sum of terms, we initialize P as
P = 0;
for k = 1:numPts
    % since the Lagrange basis polynomial is a product of terms
    % we initialize lk as
    lk = 1;
    
    % computation of lk
    for i = 1:numPts
        if i ~= k
            numer = t - T(i);
            denom = T(k) - T(i);
            lk = lk * (numer/denom);
        end
    end

    P = P + S(k) * lk;
end

%% Item 2: Displaying the polynomial P and plot of results

% polynomial P
P = simplify(P);
disp('The AIP relative to the data points is P(t) = ')
pretty(P)  % displays P decently

disp('or P(t) = ')
disp(vpa(expand(P), 4))  %  displays P having coeff in sci notation w 4 sfs 


% plot of 9 data points and interpolant P
figure
fplot(P, interval, 'k') % creates a 2-D plot of polynomial P over interval [1, 246]
hold on  % to include scatter plot with graph of P
scatter(T, S, 'filled', 'r')  % scatter plot of data points
grid  % adds grid to figure

% specify axis limits
xlim([-1 250])
ylim([160 235])

% plot labels
title('Plot of data points and interpolant $P$', ...
       'Interpreter', 'Latex', ...
       'FontSize', 16)
lgnd = legend('Interpolant $P$', 'Data Points');
xl = xlabel('day number $(t)$');
yl = ylabel('closing price of ICT stock (in pesos)');

for label = [lgnd xl yl]
    set(label,'Interpreter', 'Latex', 'Fontsize', 14)
end

%% Item 3: Estimate on a Date

t_may17 = 137;  % time variable t at May 17
price_may17 = subs(P, t_may17);  % substitute t_may17 to P
price_may17 = double(price_may17);
price_may17 = round(price_may17, 2);  % rounds price to 2 decimal places

disp(' ')
disp('the estimated closing price (in pesos) of ICT stock on May 17 is')
disp(price_may17)
disp(' ')

% computing for relative error of estimate
rel = abs(203 - price_may17) / abs(203);

disp('Hence, the relative error of the estimate is')
disp(round(rel, 4, 'significant'))  % displays rel error w 4 sf
disp('which shows that the estimate is relatively accurate.')

%% Item 4: Historical Data Points and AIP

% the historical data is imported from an excel file
% and was defined as the variable 'data'

data = table2array(data);  % to convert imported table into a matrix
T1 = data(:,1);  % time variable of historical data
S1 = data(:,2);  % stock prices of historical data

% plot of historical data points and interpolant P
figure
fplot(P, interval, 'k')
hold on
scatter(T1, S1, 'filled', 'b')
grid

% specify limits
xlim([-3 260])
ylim([165 240])

% plot labels
title('Plot of historical data points and interpolant $P$', ...
       'Interpreter', 'Latex', ...
       'FontSize', 16)
lgnd = legend('Interpolant $P$', 'Historical Data Points');
xl = xlabel('day number $(t)$');
yl = ylabel('closing price of ICT stock (in pesos)');

for label = [lgnd xl yl]
    set(label,'Interpreter', 'Latex', 'Fontsize', 14)
end

%% Item 6: Constructing new AIP (Q)

% following similar algorithm in item 1
T2 = [3 14 33 46 60 74 91 103 122 136 152 166 182 196 213 227 244 258];
S2 = [195 198.5 200.2 216 218 220 226.6 220.4 215.6 207.6 210.2 190 191 179.4 183.9 184.6 182.9 185];

numPts1 = size(T2, 2);
interval1 = [1 260];

Q = 0;
for k1 = 1:numPts1
    lk1 = 1;
    for i1 = 1:numPts1
        if i1 ~= k1
            numer1 = t - T2(i1);
            denom1 = T2(k1) - T2(i1);
            lk1 = lk1 * (numer1/denom1);
        end
    end
    Q = Q + S2(k1) * lk1;
end

% plotting new AIP (Q)
Q = simplify(Q);
disp('The AIP relative to the 18 data points is Q(t) = ')
pretty(Q)

disp('or Q(t) = ')
disp(vpa(expand(Q), 4))

% estimated value of Q at May 17
newprice_may17 = subs(Q, t_may17);
newprice_may17 = double(newprice_may17);
newprice_may17 = round(newprice_may17, 2);

disp(' ')
disp('The estimated closing price (in pesos) of ICT stock on May 17 is')
disp(newprice_may17)
disp(' ')

% computing for relative error of estimate
rel1 = abs(203 - newprice_may17) / abs(203);

disp('Thus, the relative error of the estimate is')
disp(round(rel1, 4, 'significant'))
disp('which shows that the estimate is more accurate')
disp('than the estimate obtained using the polynomial P.')

%% Item 7

% plot of polynomial Q and historical data
figure
fplot(Q, interval1, 'k')
hold on
scatter(T1, S1, 'filled', 'b')
grid

% specify limits
xlim([-3 260])
ylim([-2000 5500])

% plot labels
title('Plot of historical data points and interpolant $Q$', ...
       'Interpreter', 'Latex', ...
       'FontSize', 16)
lgnd = legend('Interpolant $Q$', 'Historical Data Points');
xl = xlabel('day number $(t)$');
yl = ylabel('closing price of ICT stock (in pesos)');

for label = [lgnd xl yl]
    set(label,'Interpreter', 'Latex', 'Fontsize', 14)
end

%% End of Program
