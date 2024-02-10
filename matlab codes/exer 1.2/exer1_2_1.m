clear
close all

syms x  % Declares x as a symbolic variable
%%
X = linspace(0, 0.5, 8);  % generate 8 equally-spaced points from 0 to 0.5
F = exp(X.^2);  % obtains f(x) values for x's in X
f = exp(x^2);  % symbolic function f(x) = e^(x^2)

numPts = length(X);
interval = [0, 0.5];

%% Calculation of functions to maximize

% f^(8)
f8Prime = diff(f, 8);  % obtains 8th order derivative of f
f8Prime = simplify(f8Prime);
disp('f^(8) =')
disp(vpa(expand(f8Prime), 4))

% plotting f8Prime
fplot(f8Prime, interval, 'r')
grid
xlim([-0.05 0.55])
ylim([1500 7800])
xl = xlabel('$x$');
yl = ylabel('$f^{(8)}(x)$');
title('Plot of $f^{(8)}(x)$ over the $x$-interval $[0, \frac{1}{2}]$', ...
    'Interpreter', 'Latex', 'Fontsize', 16)
for j = [xl yl]
    set(j, 'Interpreter', 'Latex', ...
        'Fontsize', 14)
end

% omega
omega = 1;
for i = 1:numPts
    omega = omega * (x - X(i));
end
omega = simplify(omega);
disp('omega(x) =')
disp(vpa(expand(omega), 4))

% integrating omega
integral = int(abs(omega), 0, 0.5);

%% Calculation of Error Bound
maxf8Prime = double(subs(f8Prime, 0.5));
Ebound = (maxf8Prime/factorial(numPts)) * integral;
disp('E <=')
disp(vpa(expand(Ebound), 4))

%% Calculation of AIP P
P = 0;
for k = 1:numPts
    lk = 1;
    for i = 1:numPts
        if i ~= k
            lk = lk * ((x - X(i))/(X(k) - X(i)));
        end
    end
    P = P + F(k) * lk;
end

P = simplify(P);

%% Calculation of Approximate Relative Error
Erel = Ebound/int(P, 0, 0.5);
disp('E_rel is approximately')
disp(vpa(expand(Erel),4 ))