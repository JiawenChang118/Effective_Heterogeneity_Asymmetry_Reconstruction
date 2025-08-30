function [w_IE,S_I_star] = FIC(G,w_EE,w_EI)
load('SC.mat');
SC = fln;
SC = SC./max(SC,[],'all');
N = length(SC);
% tau_E = 0.1;
tau_I = 0.01;
a_I = 615;
b_I = 177;
d_I = 0.087;
w_I = 0.7; % Scale I_0
J = 0.15; %J_NMDA
I_0 = 0.382;
% H_E_FIC = 3;%Hz
S_E_FIC = 0.164757; %Gating Varible
x_E_FIC = 0.37738; %nA

% H_E = @(x)dMFM_H(x,a_E,b_E,d_E); % f-I Curve for Excitatory Population
H_I = @(x)dMFM_H(x,a_I,b_I,d_I); % f-I Curve for Inhibitory Population

% x_I_fun = @(x,w_EI_index) -1/tau_I*(w_EI_index.*S_E_FIC-x+w_I*I_0)+H_I(x);
% [x_I_star,S_I_star] = deal(zeros(N,1));

x_I_init = 0.5;
opts = optimoptions('fsolve', ...
    'Display', 'off', 'MaxIterations', 500, ...
    'MaxFunctionEvaluations', 1000);
    x_I_star = zeros(N,1);
for i = 1:N
    x_I_fun_vec = @(x) -1/tau_I*(w_EI(i)*S_E_FIC - x + w_I*I_0) + H_I(x);
    x_I_star(i) = fsolve(x_I_fun_vec, x_I_init, opts);
end

S_I_star = tau_I * H_I(x_I_star);

w_IE = -1./S_I_star.*(x_E_FIC-w_EE.*S_E_FIC-G.*J*SC*S_E_FIC*ones(N,1)-I_0);
if min(w_IE)<0
    error('w_IE not converge.')
end
end