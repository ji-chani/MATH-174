clear
close all

syms x

%% 

numPts = 9;
% Calculation of Chebyshev polynomial of degree 9
Tnm1 = 1;  % T_0
Tn = x;  % T_1
for j = 2:numPts
    Tnp1 = 2*x*Tn - Tnm1;  % recursion formula for T_{n+1}
    Tnm1 = Tn;  % updates new T_{n-1}
    Tn = Tnp1;  % updates new T_n
end
disp('T_9 = ')
disp(vpa(expand(Tnp1), 4))

% compute for monic Chebyshev polynomial of degree 9
monicT = (2^(1-numPts)) * Tnp1;
disp('That_9 = ')
disp(vpa(expand(monicT), 4))