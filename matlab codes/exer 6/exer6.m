clc
clear
close all

syms k h x f 
%% Initializing constants
a = 0;
b = 1;
e = 10^(-8);  % convergence tolerance
RCE = e + 1; % relative computable error estimate

k(1) = 1;
k(2) = 2;
for i = 1:length(k)
    h(i) = (b-a)/2^(k(i)-1);
end

for i = 
while RCE > e


    for i = 1:2^k(i)
        x(k(i),i) = 
end



