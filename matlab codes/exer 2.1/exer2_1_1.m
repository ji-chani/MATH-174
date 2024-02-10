clear
close all

syms x
%%  Frequently used variables
Xq = linspace(-1, 1, 9);
Y = atan(x.^3) + exp(x);
numPts = length(Xq);
interval = [-1 1];

%% Calculation of AIP Q
Yq = atan(Xq.^3) + exp(Xq);
Q = 0;
for k = 1:numPts
    lk = 1;
    for i = 1:numPts
        if i ~= k
            lk = lk * ((x-Xq(i))/(Xq(k)-Xq(i)));
        end
    end
    Q = Q + Yq(k) * lk;
end
Q = simplify(Q);

%% Calculation of AIP P

% obtaining roots of T^hat_9
for j = 1:numPts
    root = cos((2*j-1)/(2*numPts) * pi);
    Xp(j) = root;
end
disp(Xp');

% final calculation for P
Yp = atan(Xp.^3) + exp(Xp);
P = 0;
for k = 1:numPts
    lk = 1;
    for i = 1:numPts
        if i ~= k
            lk = lk * ((x-Xp(i))/(Xp(k)-Xp(i)));
        end
    end
    P = P + Yp(k) * lk;
end
P = simplify(P);

%% Displaying polynomials
disp('P8(x) = ')
disp(vpa(expand(P), 4));
disp('Q8(x) = ')
disp(vpa(expand(Q), 4));

%% Plotting P8 and data points
figure
fplot(P, interval, 'k')
hold on
scatter(Xp, Yp, 'filled', 'r')
xlim([-1.2 1.2])
ylim([-0.55 3.575])
grid
lgnd = legend('$P_8$', 'Data Points');
xlp = xlabel('$x$');
ylp = ylabel('Interpolant $P_8(x)$');
for label = [lgnd xlp ylp]
    set(label, 'Interpreter', 'Latex', 'Fontsize', 14)
end
hold off

%% Plotting Q8 and data points
figure
fplot(Q, interval, 'k')
hold on
scatter(Xq, Yq, 'filled', 'b')
xlim([-1.2 1.2])
ylim([-0.55 3.55])
grid
lgnd = legend('Interpolant $Q_8$', 'Data Points');
xlp = xlabel('$x$');
ylp = ylabel('$Q_8(x)$');
for label = [lgnd xlp ylp]
    set(label, 'Interpreter', 'Latex', 'Fontsize', 14)
end

%% Plot of P8, Q8 and Y
figure
fplot(Y, interval, 'k') 
hold on
fplot(P, interval, 'r')
hold on
fplot(Q, interval, 'b')
xlim([-1.2 1.2])
ylim([-0.55 3.55])
grid

legend('$Y(x) = \arctan{x^3} + e^x$', '$P_8(x)$', '$Q_8(x)$', ...
    'Interpreter', 'Latex', 'Fontsize', 14)
xlabel('$x$', 'Interpreter', 'Latex', 'Fontsize', 14)
ylabel('function values','Interpreter', 'Latex', 'Fontsize', 14)

%% Plotting Pe and Qe
Pe = abs(P - Y);
Qe = abs(Q - Y);

figure
fplot(Pe, interval, 'b') 
hold on
fplot(Qe, interval, 'r')
grid
xlim([-1.1 1.1])
ylim([-0.1 3.75]*10^(-3))
grid

lgnd = legend('$P_e = |P_8(x) - Y(x)|$', '$Q_e = |Q_8(x) - Y(x)|$');
xlp = xlabel('$x$');
ylp = ylabel('absolute error');
for label = [lgnd xlp ylp]
    set(label, 'Interpreter', 'Latex', 'Fontsize', 14)
end

%% Error Analysis 
% due to the incapability of MATLAB to plot a directly differentiated
% f9Prime, the f9Prime is obtained using an online derivative calculator
% https://www.symbolab.com/solver/derivative-calculator/
% f9Prime is graphed using Desmos 
% from then we obtained that f9Prime is maximum whenever x=+-0.866

f9Prime = diff(Y, 9);
maxf9 = abs(double(subs(f9Prime, 0.866)));
maxY = double(subs(Y, 1));

% relative L-infty error for P8
omegaP = 1;
for i = 1:numPts
    omegaP = omegaP * (x - Xp(i));
end
omegaPprime = diff(omegaP);
CN = roots(sym2poly(omegaPprime));
CN = CN(CN >= min(interval) & CN <= max(interval));
CN = [CN' interval];
Pv = abs(double(subs(omegaP, CN)));
maxPv = max(Pv);

relP = (maxf9 * maxPv) / (factorial(9) * maxY);
disp('An upper bound for the relative L-infty error for P8 is')
disp(relP)

% relative L-infty error for Q8
omegaQ = 1;
for i = 1:numPts
    omegaQ = omegaQ * (x - Xq(i));
end
omegaQprime = diff(omegaQ);
CN = roots(sym2poly(omegaQprime));
CN = CN(CN >= min(interval) & CN <= max(interval));
CN = [CN' interval];
Qv = abs(double(subs(omegaQ, CN)));
maxQv = max(Qv);

relQ = (maxf9 * maxQv) / (factorial(9) * maxY);
disp('An upper bound for the relative L-infty error for Q8 is')
disp(relQ)

%% Calculation of AIP R8

% obtaining translated roots/nodes  [x_i = c*root + d]
intR = [0 5];  % interval for R
for j = 1:numPts
    c = (max(intR) - min(intR)) / 2;
    d = (min(intR) + max(intR)) /2 ;
    root = c * (cos((2*j-1)/(2*numPts) * pi)) + d;
    Xr(j) = root;
end

% final calculation for R
Yr = atan(Xr.^3) + exp(Xr);
R = 0;
for k = 1:numPts
    lk = 1;
    for i = 1:numPts
        if i ~= k
            lk = lk * ((x-Xr(i))/(Xr(k)-Xr(i)));
        end
    end
    R = R + Yr(k) * lk;
end
R = simplify(R);
disp('R_8(x) = ')
disp(vpa(expand(R), 4));

%% Plot for R8
figure
fplot(R, intR, 'g')
hold on
scatter(Xr, Yr, 'filled', 'k')
xlim([-0.2 5.2])
ylim([-0.5 155])
grid
lgnd = legend('Interpolant $R_8$', 'Data Points');
xlp = xlabel('$x$');
ylp = ylabel('$R_8(x)$');
for label = [lgnd xlp ylp]
    set(label, 'Interpreter', 'Latex', 'Fontsize', 14)
end

%% Calculation of L-infty norm of omega for R8
omega = 1;
for i = 1:numPts
    omega = omega * (x - Xr(i));
end

omegaPrime = diff(omega);
CN = roots(sym2poly(omegaPrime));
CN = CN(CN >= min(intR) & CN <= max(intR));
CN = [CN' intR];
Rvalues = abs(double(subs(omega, CN)));
maxRvalue = max(Rvalues);
disp('L-infty norm of omega for R8 is')
disp(maxRvalue)

%% Calculation of L-infty norm of omega for P8
omega = 1;
for i = 1:numPts
    omega = omega * (x - Xp(i));
end

omegaPrime = diff(omega);
CN = roots(sym2poly(omegaPrime));
CN = CN(CN >= min(interval) & CN <= max(interval));
CN = [CN' interval];
Pvalues = abs(double(subs(omega, CN)));
maxPvalue = max(Pvalues);
disp('L-infty norm of omega for P8 is')
disp(maxPvalue)
