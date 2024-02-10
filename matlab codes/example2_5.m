% This code is for Example 2.5 (MATH 174 - Module 2)
% The Pricing Problem

clear
close all

syms x  % declares t as a symbolic variable
%%
price = [70 82.5 137.5 187.50];
units = [1250 750 550 272.50];

numPts = length(price);
interval = [min(price) max(price)];
%% Calculate the AIP r(x)
r = 0;
for k = 1:numPts
    lk = 1;
    for i = 1:numPts
        if i ~= k
            numer = x - price(i);
            denom = price(k) - price(i);
            lk = lk * (numer/denom);
        end
    end
    r = r + units(k) * lk;
end

r = simplify(r);

% display polynomial
disp('The AIP for the given data points is r(x) =')
pretty(r)

disp('or r(x) =')
disp(vpa(expand(r), 4))


%% Graph of r
figure
fplot(r, interval, 'k');
hold on
scatter(price, units, 'filled', 'r')
grid

xlim([60 200])
ylim([200 1300])

%% Maximizing revenue R(x) = x * r(x)
R = x * r;

% differentiating R to find critical points
Rprime = diff(R);  % result is still symbolic
Rprime = sym2poly(Rprime);  % to get coefficients of Rprime

CN = roots(Rprime);  % critical points (price)
CN = CN(CN >= min(price) & CN <= max(price));  % filters critical points that are within interval
CN = [CN' [min(price) max(price)]];  % concatenate endpoints of interval (By EVT)
R_cn = double(subs(R, CN));
R_cn_max = max(R_cn);

for i = 1:length(R_cn)
    if R_cn(i) == R_cn_max
        max_index = i;
    end
end

optimal_price = CN(max_index);
disp('The optimal selling price that will result')
disp('with respect to gross revenue is')
disp(optimal_price)

