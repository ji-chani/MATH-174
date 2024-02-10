clear
close all

syms x L  % defines x and L as symbolic variables
%% Legendre roots for an AIP of deg 10
LegendreDeg = 11;  % 11 Legendre nodes to be used

L0 = 1;  % Legendre Polynomial when n=0
L(1) = x;  % to input leg poly P_1 = x to vector form L
L(2) = 3/2*x^2 - 1/2;

for n = 3:LegendreDeg
    % Recursion formula
    % Pn = (2n-1)/n *x*Pn-1 - (n-1)/n * Pn-2
    left = ((2*n-1)/n) * x * L(n-1);
    right = ((n-1)/n) * L(n-2);
    L(n) = left - right;
end

L = L(LegendreDeg);  % sets L as the degree 9 Legendre Polynomial
coeffs = sym2poly(L);
T = roots(coeffs);

numPts = length(T); % counts number of roots of Legendre Polynomials

%% Translating interval for Nodes/Abscissas
TranslateRoots = 1;  % 1 if translate, 0 if not
newInt = [-5, 3];
a = min(newInt);
b = max(newInt);

if TranslateRoots == 1
    T = ((b-a)/2) * T + mean(newInt);  % mean(newInt) = (a+b)/2
end

interpInt = [min(T) max(T)];
%% Computing AIP P
D = atan(T) - exp(T + 0.5) + T + 2;
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

disp(vpa(expand(P), 4));

%% Computation of L2 Norm
omega = 1;
for j = 1:numPts
    omega = omega * (x - T(i));
end

L2omega = int(omega^2, min(interpInt), max(interpInt));
L2omega = double(sqrt(L2omega));

disp(L2omega)

%% Absolute Error
R = abs(P - (atan(x) - exp(x+0.5) + x + 2));
figure
fplot(R, interpInt)
title('Plot of the Absolute Error')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')