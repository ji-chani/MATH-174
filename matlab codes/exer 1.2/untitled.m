clear
close all

syms x
%% Abscissas
X = linspace(0, 0.5, 8);
f = exp(x^2);

points = length(X);
interval = [0 1/2];
%% Calculation of Error Bound for E_7
% E_7 = \frac{f^{(8)}(\xi)}{8!} * \omega(x) where
% \omega(x) = \prod_{i=1}^{8}(x-x_i)

% solving for f^(8)
f8Prime = diff(f, length(X));
disp(vpa(expand(f8Prime), 4))

% solving for \omega(x)
omega = 1;
for i = X
    omega = omega * (x - i);
end

% solving for error bound
% since f^(8) is an increasing func, max(f^{(8)}) over the interval [0, 1/2]
% is f^{(8)}(1/2)
max_f8Prime = double(subs(f8Prime, 1/2));
E7_bound = (max_f8Prime/factorial(length(X))) * omega;

%% solving for bound for Error E
E_bound = int(E7_bound, 0, 1/2);
E_bound = double(E_bound);