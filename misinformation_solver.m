%% Solve the misinformation model
% define parameters


% % OUT BREAK PARAMETERS!!!!!!!!
% params.mu = 1.5; % mu > xi
% params.xi = 0.5;
% params.gamma = 0.8; % gamma > phi_N, lambda_N
% params.beta = 0.8; % beta > phi_N, lambda_N
% params.alpha = 0.2;
% params.delta = 0.2;
% params.phi_F = 0.6; % phi_F > phi_B
% params.phi_B = 0.1; % (small)
% params.phi_N = 0.1;
% params.lambda_F = 0.2; % lambda_F < lambda_B (small)
% params.lambda_B = 0.1;
% params.lambda_N = 0.1;
% params.c1 = 1.5; % initial model has c1,c2,c = 0
% params.c2 = 1.5;
% params.c = params.c1 + params.c2;


% No Outbreak
params.mu = 0.9; % mu > xi
params.xi = 0.5;
params.gamma = 0.8; % gamma > phi_N, lambda_N
params.beta = 0.8; % beta > phi_N, lambda_N
params.alpha = 0.2;
params.delta = 0.2;
params.phi_F = 0.2; % phi_F > phi_B
params.phi_B = 0.1; % (small)
params.phi_N = 0.1;
params.lambda_F = 0.6; % lambda_F < lambda_B (small)
params.lambda_B = 0.1;
params.lambda_N = 0.1;
params.c1 = 10; %11,0      initial model has c1,c2,c = 0,  for F<0.05, c1,c2=5
params.c2 = 4.5; %0,4.5
params.c = params.c1 + params.c2;


% initial conditions
% IC = [0.75; 0; 0; 0.25; 0; 0];

% outbreak onset initial conditions
IC = [0.007346; 0.135457;  0.027891; 0.434738; 0.095724; 0.298844]; % [S, E_F, E_B, F, B, N]

% time span
tspan = [0 50];
options = odeset('RelTol',1e-6,'AbsTol',1e-8,'NonNegative',1:6);
% solve
[t,Y] = ode23s(@(t,Y) misinformation_system(t,Y,params), tspan, IC, options);
% sanity check
dY0 = misinformation_system(0, [1;0;0;0;0;0], params);
disp(dY0)    % should show all zeros
%plot results
figure;
plot(t, Y, 'LineWidth', 2);
legend('S','E_F','E_B','F','B','N');
xlabel('Time');
ylabel('Population');

% plot each variable separately (same scale)
figure;
subplot(3,2,1); plot(t,Y(:,1), 'LineWidth', 2); title('S'); ylim([-0.1 max(1.1,max(Y(:,1)))]);
subplot(3,2,2); plot(t,Y(:,2), 'LineWidth', 2); title('E_F');
subplot(3,2,3); plot(t,Y(:,3), 'LineWidth', 2); title('E_B');
subplot(3,2,4); plot(t,Y(:,4), 'LineWidth', 2); title('F');
subplot(3,2,5); plot(t,Y(:,5), 'LineWidth', 2); title('B');
subplot(3,2,6); plot(t,Y(:,6), 'LineWidth', 2); title('N');
% plot on log scale
%figure;
%semilogy(t, max(abs(Y),1e-12)); legend('S','E_F','E_B','F','B','N'); xlabel('Time'); ylabel('log(|populations|)')

