%% Generate data for plot
mu    = 0.9;  % transmission rates
xi    = 0.5;
gamma = 0.8;
beta  = 0.8;
alpha = 0.2;
delta = 0.2;
phi_F = 0.6;
phi_B = 0.1;
phi_N = 0.1;
lambda_F = 0.2;
lambda_B = 0.3;
lambda_N = 0.3;

% Change to our models
c1 = linspace(0,10,100);   % you will likely use different parameter ranges
c2 = linspace(0,10,100);    % you will likely use different parameter ranges


EP = zeros(length(c1),length(c2)); % Stable EP -- one value for each choice of parameters c1 and c2

% Idea: solve the model with low ICs to get refuge EP when it exists
% (always)
for i = 1:length(c1)
    for j = 1:length(c2)
        f = @(t,u) [-mu*u(1)*u(4) - xi*u(1)*u(5);
            mu*u(1)*u(4) + gamma*u(6)*u(4) - (phi_F + phi_B + phi_N)*u(2);
            xi*u(1)*u(5) + beta*u(6)*u(5) - (lambda_F + lambda_B + lambda_N)*u(3);
            phi_F*u(2) + lambda_F*u(3) - alpha*u(4) - (c1(i)+c2(j))*(1.5*exp(-(0.25*t-1.5)^2)+0.15)*u(4);
            phi_B*u(2) + lambda_B*u(3) - delta*u(5) + c1(i)*(1.5*exp(-(0.25*t-1.5)^2)+0.15)*u(4);
            phi_N*u(2) + lambda_N*u(3) - gamma*u(6)*u(4) - beta*u(6)*u(5) + alpha*u(4) + delta*u(5) + c2(j)(1.5*exp(-(0.25*t-1.5)^2)+0.15)*u(4)];  % Define the system of ODEs with your choice of parameters from the loop
        [t,y] = ode23s(f,[0,1000],[0.007346; 0.135457;  0.027891; 0.434738; 0.095724; 0.298844]);  % solve the differential equations
        EP(i,j) = y(end,4);   % store the result-- assume that it's converged to an EP
    end
end

%% Plots

figure
hold on
surf(c1,c2,EP');    % need to transpose to get match for surf plot
colormap(flipud(parula));   % reversed parula
colorbar
xlabel('c1');
ylabel('c2');
zlabel('F');
shading interp
set(gca,'FontSize',25)

% Set 3D view
view(3);
set(gca,'FontSize',16)

