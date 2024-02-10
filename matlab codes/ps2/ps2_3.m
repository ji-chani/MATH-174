clc
clear
close all
format long

%% Initializing Constants

% integrating bounds
a = [0 1];  % lower bounds
b = [2 2];  % upper bounds

% volume of omega
V = abs(b(1)-a(1)) * abs(b(2)-a(2));

% no. of data points
N = 500;

% no. of trials
trials = 10;
t = 1;  % counter

% function and estimates
func = @(x,y) 2*y.^2 .* sin(x.*y);  % vectorized func
est = zeros(trials, 1);

while t <= trials
    % random points between bounds
    X = (b(1)-a(1))*rand(N,1) + a(1);
    Y = (b(2)-a(2))*rand(N,1) + a(2);

    % function values
    F = 2*Y.^2 .* sin(X.*Y);

    % t^th trial estimate
    est(t) = V * mean(F);
    
    %% Plotting
    if t == 1
        [X_3d, Y_3d] = meshgrid(a(1):0.1:b(1), a(2):0.1:b(2));
        Z_3d = 2*Y_3d.^2 .* sin(X_3d.*Y_3d);  % plot of f(x,y) in 3d

        surf(X_3d, Y_3d, Z_3d)
        colormap("pink")
        hold on
        scatter3(X,Y,F,'ko','filled')
        
        % labels
        legend('$f(x,y) = 2y^2 \sin{(xy)}$', 'sample points', ...
            'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$y$', 'Interpreter', 'latex')
        zlabel('$z$', 'Interpreter', 'latex')
        zlim([-8 10])
        hold off
    end
    t = t+1;
end
disp('Estimates = ')
disp(est)

%% Relative Error
exact = integral2(func, a(1), b(1), a(2), b(2));  % int(func, xmin, xmax, ymin, ymax)

RelError = abs((exact-est)/exact) * 100;
disp('Relative errors (in %) = ')
disp(RelError)

%% Average estimate and Relative Error
est_ave = mean(est);
relE = abs((exact-est_ave)/exact) * 100;
disp('average estimate =')
disp(est_ave)
disp('with relative error (in %) =')
disp(relE)