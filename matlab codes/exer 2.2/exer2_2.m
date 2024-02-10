clear
close all

syms x L  % Declares x and L as symbolic variables (L is a symbolic array)
%% Obtaining Chebyshev and Legendre Roots

% Chebyshev
ChebDeg = 9;  % Chebyshev Degree
Tp = zeros(1, ChebDeg);  % chebyshev nodes
for j = 1:ChebDeg
    root = cos((2*j-1)/(2*ChebDeg) * pi);
    Tp(j) = root;
end

% Legendre
LegDeg = 9;  % Legendre Degree
L0 = 1;  % L_0 = 1
L(1) = x;  % L_1 = x
L(2) = 3/2*x^2 - 1/2;
for n = 3:LegDeg
    left = ((2*n-1)/n)*x*L(n-1);
    right = ((n-1)/n)*L(n-2);
    L(n) = left - right;
end

L = L(LegDeg);
Ts = roots(sym2poly(L));  % Legendre nodes

numPts = length(Ts);  % counts number of interpolatory abscissas

%% Compute for AIP P Q and S

% AIP P
Kp = atan(Tp.^3) + exp(Tp);
P = 0;
for k = 1:numPts
    lk = 1;
    for i = 1:numPts
        if i ~= k
            numer = x - Tp(i);
            denom = Tp(k) - Tp(i);
            lk = lk * (numer/denom);
        end
    end
    P = P + Kp(k) * lk;
end

% AIP Q
Tq = linspace(-1, 1, 9);  % equally spaced points
Kq = atan(Tq.^3) + exp(Tq);
Q = 0;
for k = 1:numPts
    lk = 1;
    for i = 1:numPts
        if i ~= k
            numer = x - Tq(i);
            denom = Tq(k) - Tq(i);
            lk = lk * (numer/denom);
        end
    end
    Q = Q + Kq(k) * lk;
end

% AIP S
Ks = atan(Ts.^3) + exp(Ts);
S = 0;
for k = 1:numPts
    lk = 1;
    for i = 1:numPts
        if i ~= k
            numer = x - Ts(i);
            denom = Ts(k) - Ts(i);
            lk = lk * (numer/denom);
        end
    end
    S = S + Ks(k) * lk;
end

disp('S8(x) = ')
disp(vpa(expand(S), 4))

%% Plotting
Y = atan(x^3) + exp(x);
interval = [-1, 1];

% S8, data points, and Y
figure
fplot(Y, interval, 'k')
hold on
fplot(S, interval, 'r')
hold on
scatter(Ts, Ks, 'filled', 'b')
xlim([-1.1 1.1])
ylim([-1 5.3])
grid
lgnd = legend( ...
    '$Y(x) = \arctan{(x^3)} + e^x$', ...
    'Interpolant $S_8(x)$', ...
    '"Legendre nodes"');
xlp = xlabel('$x$');
ylp = ylabel('$y$');
for label = [lgnd xlp ylp]
    set(label, 'Interpreter', 'Latex', 'Fontsize', 14)
end

% P8, Q8, S8, and K
figure
fplot(Y, interval, 'k')
hold on
fplot(P, interval, 'g')
hold on
fplot(Q, interval, 'b')
hold on
fplot(S, interval, 'r')
xlim([-1.1 1.1])
ylim([-1 5.3])
grid
lgnd = legend( ...
    '$Y(x) = \arctan{(x^3)} + e^x$', ...
    'Interpolant $P_8(x)$', ...
    'Interpolant $Q_8(x)$', ...
    'Interpolant $S_8(x)$');
xlp = xlabel('$x$');
ylp = ylabel('$y$');
for label = [lgnd xlp ylp]
    set(label, 'Interpreter', 'Latex', 'Fontsize', 14)
end

%% Plotting absolute errors
Pe = abs(P-Y);
Qe = abs(Q-Y);
Se = abs(S-Y);

figure
fplot(Pe, interval, 'g')
hold on
fplot(Qe, interval, 'b')
hold on
fplot(Se, interval, 'r')
xlim([-1.1 1.1])
ylim([-0.0001 .004])
grid
lgnd = legend('$P_e$', '$Q_e$', '$S_e$');
xlp = xlabel('$x$');
ylp = ylabel('$y$');
for label = [lgnd xlp ylp]
    set(label, 'Interpreter', 'Latex', 'Fontsize', 14)
end

%% Upper bound for Relative L2 Error

% L2 norm of omegas
omegaP = 1;
for i = 1:numPts
    omegaP = omegaP * (x - Tp(i));
end

omegaQ = 1;
for i = 1:numPts
    omegaQ = omegaQ * (x - Tq(i));
end

omegaS = 1;
for i = 1:numPts
    omegaS = omegaS * (x - Ts(i));
end

omegaPnorm = double(sqrt(int(omegaP^2, min(Tp), max(Tp))));
omegaQnorm = double(sqrt(int(omegaQ^2, min(Tq), max(Tq))));
omegaSnorm = double(sqrt(int(omegaS^2, min(Ts), max(Ts))));

% L2 Norm of Function Y(x) = atan(x^3)+exp(x)
y9Prime = diff(Y, 9);
maxy9Prime = double(abs(subs(y9Prime, 0.866)));

Ynorm = double(sqrt(int(Y^2, -1, 1)));

% Computation of Relative Errors
factor = maxy9Prime/abs(factorial(numPts));
relP = (factor * omegaPnorm) / Ynorm;
relQ = (factor * omegaQnorm) / Ynorm;
relS = (factor * omegaSnorm) / Ynorm;

disp('relp = ')
disp(vpa(relP, 4))
disp('relq = ')
disp(vpa(relQ, 4))
disp('rels = ')
disp(vpa(relS, 4))

%% Root Translation
TranslateRoots = 1;  % Indicator if I wish to scale and translate the abscissas
newInt = [0, 5];
a = newInt(1);
b = newInt(2);

if TranslateRoots == 1
    % scaled Chebyshev Nodes
    Tr = ((b-a)/2)*Tp + mean(newInt);  % mean(newInt) = (a+b)/2
    % scaled Legendre Nodes
    Tt = ((b-a)/2)*Ts + mean(newInt);
end

interpIntR = [min(Tr) max(Tr)];
interpIntT = [min(Tt) max(Tt)];

%% Computation of AIP R and T
% Translated Chebyshev
Kr = atan(Tr.^3) + exp(Tr);
R = 0;
for k = 1:numPts
    lk = 1;
    for i = 1:numPts
        if i ~= k
            numer = x - Tr(i);
            denom = Tr(k) - Tr(i);
            lk = lk * (numer/denom);
        end
    end
    R = R + Kr(k) * lk;
end

% Translated Legendre
Kt = atan(Tt.^3) + exp(Tt);
T = 0;
for k = 1:numPts
    lk = 1;
    for i = 1:numPts
        if i ~= k
            numer = x - Tt(i);
            denom = Tt(k) - Tt(i);
            lk = lk * (numer/denom);
        end
    end
    T = T + Kt(k) * lk;
end

disp('T8(x) = ')
disp(vpa(expand(T), 4))
%% Plotting AIP T, R, and K
figure
fplot(Y, newInt, 'k')
hold on
fplot(R, interpIntR, 'b')
hold on
fplot(T, interpIntT, 'r')
grid
xlim([-0.1 5.1])
ylim([0 150])
lgnd = legend('$Y(x) = \arctan{(x)} + e^x$', ...
    'Interpolant $R_8(x)$', ...
    'Interpolant $T_8(x)$');
xlp = xlabel('$x$');
ylp = ylabel('$y$');
for label = [lgnd xlp ylp]
    set(label, 'Interpreter', 'Latex', 'Fontsize', 14)
end

%% Re and Te
figure
fplot(abs(R-Y), newInt, 'b')
hold on
fplot(abs(T-Y), newInt, 'r')
grid
legend('$R_e(x)$', '$T_e(x)$', 'Interpreter', 'Latex')