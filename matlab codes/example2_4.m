%  This code is for Example 2.4 (MATH 174 - Module 2)
clear
close all

syms t  % declares t as a symbolic variable
%%  Initializes abscissas from data points
T = [0 10 25 50 100]; % day number
D = [0 8 4251 20677 357];  % infected individuals

numPts = size(T, 2);
interval = [min(T) max(T)];

%% Calculating for Lagrange form of AIP (P)
P = 0;

for k = 1:numPts
    lk = 1;

    for i = 1:numPts
        if i ~= k
            numer = t - T(i);
            denom = T(k) - T(i);
            lk = lk * (numer/denom);
        end
    end
    P = P + D(k) * lk;
end

%% Postprocessing (Display of result)
P = simplify(P);
disp('The AIP relative to the given data points is P(t) = ')
pretty(P)

disp('or P(t) = ')
disp(vpa(expand(P), 4))

%% Weekly Infections
t = 7:7:100;

% to substitute each of elements of t to polynomial P
P_t = subs(P, t);
P_t = double(P_t);

disp(P_t)

int = [min(T):1:max(T)];
mx = max(subs(P, int));  % max among function values of P(t) where t = [1, 100]
mx = double(mx);
disp('The maximum number of people infected')
disp('at any given day within the 100-day period is')
disp(mx)
    