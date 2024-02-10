clc
clear
close all
format long

syms x y
%% Initialization of constants

% integration bounds
a1 = -5;
b1 = 4;
a2 = -4;
b2 = 5;

% volume
V = abs(b1-a1) * abs(b2-a2);

% no. of data points
N1 = 10^2;
N2 = 10^6;

% no. of trials
trials = 1;
t = 1;  % counter

% function and estimate
f = y*sin(x) - x*cos(y);
est = zeros(2, trials);  % placeholder for estimates

while t <= trials
    data1 = rand(N1,2);
    data2 = rand(N2,2);
    
    % for N1
    data1(:,1) = (b1-a1)*data1(:,1) + a1;  % x1
    data1(:,2) = (b2-a2)*data1(:,2) + a2;  % y1

    % for N2
    data2(:,1) = (b1-a1)*data2(:,1) + a1;  % x2
    data2(:,2) = (b2-a2)*data2(:,2) + a2;  % y2

    % coordinates
    x1 = data1(:,1);
    y1 = data1(:,2);
    x2 = data2(:,1);
    y2 = data2(:,2);

    % function values
    F1 = y1.*sin(x1) - x1.*cos(y1);
    F2 = y2.*sin(x2) - x2.*cos(y2);

    est(1,t) = V * mean(F1);
    est(2,t) = V * mean(F2);

    t = t+1;
end

%% Plotting Sample Points
[X, Y] = meshgrid(-5:0.5:4, -4:0.5:5);
Z = Y.*sin(X) - X.*cos(Y);
z = zeros(N1, 1);
for i = 1:N1
    z(i) = double(subs(f, [x, y], [x1(i), y1(i)]));
end

tiledlayout(1,2)
nexttile
surf(X, Y, Z)
colormap("pink")
hold on
scatter3(x1, y1, z, 'k', 'filled')
hold off

nexttile
surf(X, Y, Z)
colormap("pink")
hold on
scatter3(x1, y1, z, 'k', 'filled')
hold off
shading interp
view(2)