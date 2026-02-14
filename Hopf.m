% Parameters
clear
r_values = linspace(0.8, 1.2, 100);  % Range of r values around bifurcation
tspan = [0 1000];                      % Time span for integration
num_points = 500;                     % Number of points to capture at the end

% Initial conditions
x0 = 0.75;  % Initial x value
y0 = 0.75;  % Initial y value

options = odeset('RelTol',1e-8, 'AbsTol',1e-10);  % Adjust tolerances as needed

% Prepare figure
figure;
hold on;
xlabel('r');
ylabel('x');
zlabel('y');
title('Hopf Bifurcation');
grid on;

% Loop over r values to observe steady-states and limit cycles
for i = 1:length(r_values)
    r = r_values(i);
    
    % Define the system of differential equations for the current r
    dxdt = @(t, y) [r - y(1)*y(2)^2; y(1)*y(2)^2 - y(2)];
    
    % Solve the ODE
    [~, Y] = ode45(dxdt, tspan, [x0, y0]);
    
    % Select the last `num_points` values to check for limit cycles
    Y_tail = Y(end-num_points+1:end, :);
    
    % Plot the selected points with the same r value
    plot3(r * ones(num_points, 1), Y_tail(:, 1), Y_tail(:, 2), 'b.');
end

% Set 3D view
view(3);
ylim([0,3])
zlim([0,3])
set(gca,'FontSize',16)

