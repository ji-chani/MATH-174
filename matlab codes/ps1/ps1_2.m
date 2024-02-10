clc
clear
close all

syms x

%% Function and interval
f = sin(x^2);
interval = [-1, 1];

%% L-infty norm of AIP using monic Chebyshev polynomial roots
ChebDeg = 10;
omega_norm = 2^(1-ChebDeg);
f10Prime = diff(f, 10);

% Chebyshev nodes (not important in calculation)
X = zeros(1, ChebDeg);
for i = 1:ChebDeg
    root = cos((2*i-1)/(2*ChebDeg) * pi);
    X(i) = root;
end
Y = sin(X.^2);

% max value of f10Prime occurs at x=+-1
maxf10Prime = double(abs(subs(f10Prime, 1)));

% L-infty error of AIP
AIP_norm = (maxf10Prime/factorial(ChebDeg)) * omega_norm;
disp(AIP_norm)

%% L-infty norm of Piecewise Linear Interpolant
% 1/8 h^2 max(f2Prime)
f2Prime = diff(f, 2);

% max value of abs(f2Prime) occurs at x=+-1
maxf2Prime = double(abs(subs(f2Prime, 1)));

%% Plots of f2Prime and f10Prime together
t = tiledlayout(2, 1);
t.TileSpacing = 'compact';

% Tile 1
nexttile
fplot(abs(f2Prime), interval, 'k')
l1 = legend("$|f''(x)|$");
x1 = xlabel('$x$');
y1 = ylabel('$y$');

% Tile 2
nexttile
fplot(abs(f10Prime), interval, 'r')
l2 = legend('$|f^{(10)}(x)|$');
x2 = xlabel('$x$');
y2 = ylabel('$y$');
% labels
for i = [l1 l2 x1 x2 y1 y2]
    set(i, 'Interpreter', 'Latex', 'Fontsize', 14)
end

%% solving for increment h
% 1/8 h^2 max(f2Prime) = AIP_norm
h = abs(sqrt((8 * AIP_norm)/maxf2Prime));
disp(h)