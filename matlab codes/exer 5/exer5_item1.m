clc
clear
format long

syms x
%% Initializing Necessary Variables
a = 2;
N = 10^6;
f = 1/(log(x));  % f(x) = 1/ln(x) -- log in MATLAB is ln


nP4 = 4;  % num of points for 4-point CNCQR and ONCQR (n+1)
nP5 = 5;  % num of points for 5-point CNCQR and ONCQR (n+1)

%% Solving for Increments, Data Points, and Function Values

% increments h
h = zeros(4, 1);  % increments for 4pt and 5-point CNCQR,ONCQR
h(1) = (N-a)/(nP4-1);  % h = (b-a)/n where n+1=4
h(2) = (N-a)/(nP5-1);  % h = (b-a)/n where n+1=5
h(3) = (N-a)/(nP4+1);  % h = (b-a)/(n+2) where n+1=4
h(4) = (N-a)/(nP5+1);  % h = (b-a)/(n+2) where n+1=5

% data points x_i
X = zeros(3, nP5);
for i = 1:nP4
    X(1,i) = a + (i-1)*h(1);  % 4pt CNCQR
    X(3,i) = a + i*h(3);  % 4pt ONCQR
end
for i = 1:nP5
    X(2,i) = a + (i-1)*h(2);  % 5pt CNCQR
    X(4,i) = a + i*h(4);  % 5pt ONCQR
end

% solving F(X) where F(x) = 1/ln(x)
F = (log(X)).^(-1);

%% Estimates
Est = zeros(4,1); 

% 4pt closed NCQR (3/8 Rule)
Est(1) = (N-a)/8 * (F(1,1) + 3*F(1,2) + 3*F(1,3) + F(1,4));

% 5pt closed NCQR (Boole's Rule)
Est(2) = (N-a)/90 * (7*F(2,1) + 32*F(2,2) + 12*F(2,3) + 32*F(2,4) +7*F(2,5));

% 4pt open NCQR
Est(3) = (N-a)/24 * (11*F(3,1) + F(3,2) + F(3,3) + 11*F(3,4));

%% 5pt open NCQR
% integration AIP of order 4 (P)
P = 0;
for k = 1:nP5
    lk = 1;
    for i = 1:nP5
        if i ~= k
            numer = x - X(4,i);
            denom = X(4,k) - X(4,i);
            lk = lk * (numer/denom);
        end 
    end
    P = P + F(4,k) * lk;
end

Est(4) = double(int(P, a, N));

disp('The estimates are')
disp(Est)

%% Calculation of Relative Error
exact = 78498;
relError = abs((Est-exact)/exact);

disp("The relative errors are")
disp(relError * 100)