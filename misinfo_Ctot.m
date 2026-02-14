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

tolF = 0.05;
epsilon = 5e-4;

c1_values = linspace(0, 12, 50);   % high resolution for c1
valid_c1 = [];
valid_c2 = [];

tspan = [0 60];                     % shorter time is faster
IC = [0.007346; 0.135457; 0.027891; 0.434738; 0.095724; 0.298844];

options = odeset('RelTol',1e-5,'AbsTol',1e-5);   % MUCH faster

for i = 1:length(c1_values)

    params.c1 = c1_values(i);

    % --- bracket c2 by scanning a coarse grid ---
    c2_scan = linspace(0, 5, 50);
    Fvals = zeros(size(c2_scan));

    for k = 1:length(c2_scan)
        params.c2 = c2_scan(k);
        [~,Y] = ode15s(@(t,Y) misinformation_system_C(t,Y,params), tspan, IC, options);
        Fvals(k) = Y(end,4);
    end

    % Find intervals where F crosses tolF
    idx = find(diff(sign(Fvals - tolF)) ~= 0);

    if isempty(idx)
        continue
    end

    % % Now refine c2 in that specific interval using bisection
    % c2_left = c2_scan(idx);
    % c2_right = c2_scan(idx+1);
    % 
    % for iter = 1:20   % 20 bisection iterations ≈ perfect accuracy
    %     c2_mid = (c2_left + c2_right)/2;
    %     params.c2 = c2_mid;
    %     [~,Y] = ode15s(@(t,Y) misinformation_system_C(t,Y,params), tspan, IC, options);
    %     if Y(end,4) > tolF
    %         c2_left = c2_mid;
    %     else
    %         c2_right = c2_mid;
    %     end
    % end

    % c2_final = (c2_left + c2_right)/2;

    valid_c1(end+1) = params.c1;
    valid_c2(end+1) = c2_scan(idx);
end

figure; hold on;
plot(valid_c1, valid_c2, 'b-', 'LineWidth', 2)
xlabel('c_1')
ylabel('c_2')
title('Boundary curve where F = 0.05')
grid on


%%
function F_end = F_function(c1, c2, tspan, IC, params)
% F_function  Integrate the misinformation ODEs with given c1,c2 and
% return the final value of F (component 4) at tspan(end).
%
% Usage:
%  F_end = F_function(c1, c2, tspan, IC, params)
%
% Inputs:
%  c1, c2    - scalar control parameters
%  tspan     - [t0 tf] integration interval (e.g. [0 60])
%  IC        - 6x1 initial condition [S; EF; EB; F; B; N]
%  params    - struct containing model parameters:
%              params.mu, params.xi, params.gamma, params.beta,
%              params.alpha, params.delta,
%              params.phi_F, params.phi_B, params.phi_N,
%              params.lambda_F, params.lambda_B, params.lambda_N
%
% Output:
%  F_end     - scalar F(tfinal)

% choose solver options here (adjust tolerances as needed)
options = odeset('RelTol',1e-7,'AbsTol',1e-9);

% pack c1,c2 into params so the ODE can read them (optional)
params.c1 = c1;
params.c2 = c2;

% integrate (use ode15s as before)
[~, Y] = ode15s(@(t,y) misinformation_system_C(t,y,params), tspan, IC, options);

% return final F (component 4)
F_end = Y(end,4);

end

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

IC = [0.007346; 0.135457;  0.027891; 0.434738; 0.095724; 0.298844]; % [S, E_F, E_B, F, B, N]
% time span
tspan = [0 50];

% Range of c2 to search
c1_vals = linspace(0,15,2000);   % increase resolution if needed
F_vals = zeros(size(c2_vals));

% Evaluate F across this slice
for i = 1:length(c1_vals)
    F_vals(i) = F_function(c1_vals(i), 0, tspan, IC, params);  % <-- your F(c1, c2)
end

% 1. Maximum possible value of F at c1 = 0
[maxF, idxMax] = max(F_vals);
c1_at_maxF = c1_vals(idxMax);

fprintf('Max F at c2 = 0 is %.6f at c1 = %.6f\n', maxF, c1_at_maxF);

% 2. Closest value below 0.05
target = 0.05;
[~, idxClosest] = min(abs(F_vals - target));

closestF = F_vals(idxClosest);
c1_closest = c1_vals(idxClosest);

fprintf('Closest F to %.3f is %.6f at c1 = %.6f\n', target, closestF, c1_closest);
