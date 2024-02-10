clear
close all

syms x
%% Solving for omega(x)
% = (x-0)(x-pi/2)(x-pi)(x-2pi)

omega = 1;
T = [0 pi/2, pi, 2*pi];
for i = 1:length(T)
    omega = omega * (x - T(i));
end

disp(omega)

%% maximizing omega using derivatives

omegaPrime = diff(omega);  % take derivative of omega
coeff = sym2poly(omegaPrime);  % obtains coeff of omegaPrime
CN = roots(coeff);  % obtains critical points of omegaPrime

CN = CN(CN >= min(T) & CN <= max(T));  % filter entries of CN that are in [min(T), max(T)]
CN = [CN' [min(T) max(T)]];  % includes endpoints of T in CN

% substitute CN to omega
omega_cn = abs(double(subs(omega, CN)));

disp("The values of omega, respectively, at x = ")
disp(CN)
disp('is')
disp(omega_cn)

% obtaining max omega (upper bound)

% % gets index of maximum value of omega
% for i = 1:length(omega_cn)
%     if omega_cn(i) == mx
%         max_index = i;
%     end
% end

disp('Therefore, the maximum value of omega over the interval [0, 2pi] is')
disp(max(omega_cn))


