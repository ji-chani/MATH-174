clc
clear
close all
syms S
syms x real
%% Abscissas
numNodes = 7;
T = linspace(-3, 3, numNodes);
interval = [min(T) max(T)];
F = zeros(1,numNodes);  % F(x)
f = 1/(sqrt(2*pi)) * exp(-(x^2/2));  % symbolic f(x)
for i = 1:numNodes
    den = sqrt(2*pi);
    n = (T(i)^2/2);
    F(i) = (1/den)*exp(-n);
end

disp(F')

%% Calculation of linear interpolants S_i
for i = 1:numNodes-1
    a = F(i);
    b = (F(i+1)-F(i))/(T(i+1)-T(i));
    S(i) = a + b*(x-T(i));
end
S = S';
disp(vpa(expand(S), 4))
%% Plotting the Piecewise Linear Interpolant s
s = piecewise( ...
    (T(1)<=x)&(x<=T(2)), S(1), ...
    (T(2)<=x)&(x<=T(3)), S(2), ...
    (T(3)<=x)&(x<=T(4)), S(3), ...
    (T(4)<=x)&(x<=T(5)), S(4), ...
    (T(5)<=x)&(x<=T(6)), S(5), ...
    (T(6)<=x)&(x<=T(7)), S(6));

figure
hold on
fplot(f, interval, 'r')
fplot(s, 'b')
scatter(T, F, 'filled', 'k')
grid
xlim([-3.1 3.1])
ylim([-0.1 0.5])
lgnd = legend('$f(x)=\frac{1}{\sqrt{2\pi}}e^{-\frac{x^2}{2}}$', ...
    'Interpolant $s(x)$', ...
    'Abscissas', 'Interpreter', 'latex');
xlabel('$x$','Interpreter', 'latex')
ylabel('$y$','Interpreter', 'latex')
set(lgnd, 'Fontsize', 14)

%% Calculation of s(0.5) and s(1.75)
s50 = double(subs(S(4), 0.5));
s175 = double(subs(S(5), 1.75));
disp('s(0.5) =')
disp(s50)
disp('s(1.75) = ')
disp(s175)

%% Relative error of s(0.5)
exact50 = 0.352065326764300;
exact175 = 0.086277318826512;
rel50 = abs(exact50-s50)/abs(exact50);
rel175 = abs(exact175-s175)/abs(exact175);
disp('The relative error of s(0.5) is')
disp(rel50)
disp('or approximately (in percent)')
disp(vpa(rel50*100, 4))
disp('and the relative error of s(1.75) is')
disp(rel175)
disp('or approximately (in percent)')
disp(vpa(rel175*100, 4))
%% Item 2: Finding upper bound for increment h
f2Prime = diff(f, 2);

% plot f2Prime to find max
figure
fplot(abs(f2Prime), interval, 'r', 'LineWidth', 1.5)
grid
xlim([-3.1 3.1])
ylim([-0.1 0.5])
xlabel('$x$','Interpreter', 'latex')
ylabel('$|f^{(2)}(x)|$','Interpreter', 'latex')

% calculation of upper bound
maxf2Prime = double(subs(abs(f2Prime), 0));
h_ub = abs(sqrt((8*10^(-6)/maxf2Prime)));
disp('Upper bound for increment is')
disp(h_ub)