clear
close all

syms x L
%% Legendre roots for an AIP of degree 8
LegendreDegree = 11;  % 
L0 = 1;  % so that ith entry of vector L is the degree i Legendre Polynomial
L(1) = x;
L(2) = (3/2)*x^2 - (1/2);

for i = 3:LegendreDegree
    L(i) = ((2*i-1)/i)*x*L(i-1) + ((i-1)/i)*L(i-2);  % Recursion formula
end

L = L(LegendreDegree);  % sets L as the degree 9 Legendre Polynomial
coeffs = sym2poly(L);
T = roots(coeffs);

numPts = length(T);  % counts number of roots of Legendre Polynomials

%% Root Translator
TranslateRoots = 1;  % Indicator if I wish to scale and translate the abscissas
newInt = [-5, 3];
a = newInt(1);
b = newInt(2);

if TranslateRoots == 1
    T = ((b-a)/2)*T + mean(newInt);  % mean(newInt) = (a+b)/2
end

interpInt = [min(T) max(T)];

%% Compute the AIP
D = atan(T) - exp(T+0.5) + T + 2;  % the function to approximate
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

omega = 1;
for i = 1:length(T)
    omega = omega * (x-T(i));
end

% computing L-infty norm of omega
omegaP = diff(omega);
coeffsP = sym2poly(omegaP);
CN = roots(coeffsP);
Fvals = double(subs(omega, [CN', min(interpInt), max(interpInt)]));  % [min(T) max(T)] is the min and max of interpInt

% computing for L-2 norm of omega
L2omega = int(omega^2, min(interpInt), max(interpInt));
L2omega = double(sqrt(L2omega));

%%
R = abs(P - (atan(x) - exp(x+0.5) + x + 2));

figure
fplot(R, interpInt,'k', 'Linewidth', 2)
title('Plot of the Absolute Error')
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$y$', 'Interpreter', 'Latex')
