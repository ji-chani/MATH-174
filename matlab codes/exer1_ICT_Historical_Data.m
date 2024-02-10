close all

t = data(:,1);
prices = data(:,2);

t = t';
prices = prices';

%% 

figure
fplot(P, interpInt, 'k')
hold on
scatter(t, prices, 'filled', 'r')
grid

title('Plot of data points and its interpolant')
lgnd = legend('Interpolant $P$', 'Data Points', 'Interpreter', 'Latex');
x_l = xlabel('time variable ($t$)', 'Interpreter', 'Latex');
y_l = ylabel('closing stock price (in pesos)', 'Interpreter', 'Latex');

for k = [lgnd x_l y_l]
    set(k, 'Fontsize', 14)
end