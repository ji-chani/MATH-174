clear
close all

syms x  % Declares x as a symbolic variable
%%
X = [0.00480367989919189 0.0421325969243644 0.495196320100808 0.457867403075636 0.111107441745099 0.388892558254901 0.201227419495968 0.298772580504032];
F = exp(X.^2);  % obtains f(x) values for x's in X
f = exp(x^2);  % symbolic function f(x) = e^(x^2)

numPts = length(X);
interval = [0, 0.5];

%% Calculation of functions to maximize

% f^(8)
f8Prime = diff(f, 8);  % obtains 8th order derivative of f
f8Prime = simplify(f8Prime);

% omega
omega = 1;
for i = 1:numPts
    omega = omega * (x - X(i));
end
omega = simplify(omega);

% integrating omega
integral = int(abs(omega), 0, 0.5);

%% Calculation of Error Bound
maxf8Prime = double(subs(f8Prime, 0.5));
Ebound = (maxf8Prime/factorial(numPts)) * integral;

%% Calculation of AIP P
Q = 0;
for k = 1:numPts
    lk = 1;
    for i = 1:numPts
        if i ~= k
            lk = lk * ((x - X(i))/(X(k) - X(i)));
        end
    end
    Q = Q + F(k) * lk;
end

Q = simplify(Q);

%% Calculation of Approximate Relative Error
Erel = Ebound/int(Q, 0, 0.5);
disp('E_rel is approximately')
disp(vpa(expand(Erel),4 ))