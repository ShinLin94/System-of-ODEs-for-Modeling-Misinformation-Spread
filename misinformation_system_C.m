%% Create the misinformation model system of differential equations

function dYdt = misinformation_system(t, Y, params)

    % define variables and parameters
    S = Y(1);
    E_F = Y(2);
    E_B = Y(3);
    F = Y(4);
    B = Y(5);
    N = Y(6);

    mu = params.mu;
    xi = params.xi;
    gamma = params.gamma;
    beta = params.beta;
    alpha = params.alpha;
    delta = params.delta;

    phi_F = params.phi_F;
    phi_B = params.phi_B;
    phi_N = params.phi_N;
    lambda_F = params.lambda_F;
    lambda_B = params.lambda_B;
    lambda_N = params.lambda_N;

    c = params.c;
    c1 = params.c1;
    c2 = params.c2;

    % system of differential equations
    dSdt = -mu*S*F - xi*S*B;
    dEFdt = mu*S*F + gamma*N*F - (phi_F + phi_B + phi_N)*E_F;
    dEBdt = xi*S*B + beta*N*B - (lambda_F + lambda_B + lambda_N)*E_B;
    dFdt = phi_F*E_F + lambda_F*E_B - alpha*F - c*(1.5*exp(-(0.25*t-1.5)^2)+0.15)*F;
    dBdt = phi_B*E_F + lambda_B*E_B - delta*B + c1*(1.5*exp(-(0.25*t-1.5)^2)+0.15)*F;
    dNdt = phi_N*E_F + lambda_N*E_B - gamma*N*F - beta*N*B + alpha*F + delta*B + c2*(1.5*exp(-(0.25*t-1.5)^2)+0.15)*F;

    % return the solutions as a column vector
    dYdt = [dSdt; dEFdt; dEBdt; dFdt; dBdt; dNdt];

end
