%% Parameters setting

load('SC.mat');
SC = fln;
N = length(SC);
T = 50000; %s;
dt = 0.01; %10ms;
w = linspace(0.8, 1.3, N)' * 0.5; % Local excitatory recurrent; %0.5327
I = linspace(1, 1.167, N)' * 0.3; % nA, the overall effective external input; %0.3266

G = 0.65;
sigma = 0.01; % nA, the noise amplitud.
J = 0.2609;
tau = 0.1;
gamma = 0.641;

%% Setting f-I function H and dH
a = 270; %n/C (/nC?)
b = 108; % Hz
d = 0.154; % s
H = @(x)dMFM_H(x,a,b,d); % f-I Curve
dH = @(x) - 270./(exp(2079/125 - (2079*x)./50) - 1) -...
     (2079.*exp(2079/125 - (2079*x)./50).*(270*x - 108))./(50*(exp(2079/125 - (2079*x)./50) - 1).^2);

%% Reduced Wong-Wang Model Simulation

[S,eta] = dMFM(SC, dt, T, w, I, G, sigma);
S_star = mean(S,2);
x_star = w.*J.*S_star+G.*J.*SC*S_star+I;

dH_val = dH(x_star); % reveal the value of dH/dt at x_star

%% Jacobian Matrix Calculation

Jacob = zeros(N);
y_ana_inv = zeros(N,1);
for i = 1:N
    y_ana_inv(i) = gamma*G*J*(1-S_star(i))*dH(x_star(i));
    for j = 1:N
        if i == j
            Jacob(i,j) = -1/(tau*(1-S_star(i)))+w(i)*gamma*J*(1-S_star(i))*dH(x_star(i));
        else
            Jacob(i,j) = gamma*G*J*(1-S_star(i))*SC(i,j)*dH(x_star(i));
        end
    end
end
clear i j 

%% Linear Reconstruct Jacobian Matrix with S
Jacob_est = LinearReconst(S,dt); % Estimating Jacobian Matrix

%% Reconstruct spatial properties with Jacobian Matrix
SC_sym = (SC + SC')/2;
[y_st,w_recon,C_recon] = RevealHHetero1(SC_sym,Jacob_est,S_star,tau,gamma,G, J);
C_recon = C_recon-diag(diag(C_recon));

%% Evaluation of Jacobian Matrix J and SC Estimation

[SSE_J,Corr_J,Corr_nonzero_J] = EstimationJacobianPlotting(Jacob,Jacob_est);
[SSE_SC,Corr_SC,Corr_nonzero_SC] = EstimationMatrixPlotting(SC,C_recon);

%% Evaluation of Partial Firing Rate dH/dx Estimation

dH_st = 1./(gamma*G*J.*(1-S_star).*y_st);

figure(5)
scatter(dH_val,dH_st);
set(gca,'box','off');
xlabel('Ground Truth dH/dx_i^*');
ylabel('Estimated dH/dx_i^*');
xticks([0 300]);
yticks([0 300]);
alpha(0.8);
saveas(gcf,'hPerformance.png');

%% Evaluation of Local excitatory recurrent w Estimation

figure(6)
scatter(w,w_recon);
set(gca,'box','off');
hold on
plot([0.4,max(w)],[0.4,max(w)],'--','Color','k');
hold off
xlabel('Ground Truth w_i');
ylabel('Estimated w_i');
xticks([0.4 0.65]);
yticks([0.4 0.65]);
alpha(0.8);
saveas(gcf,'wPerformance.png');

%% Evaluation of Input Current x* Estimation

x_st = zeros(size(dH_st));
for i = 1:length(dH_st)
    x_initial_guess = 0.5;
    x_solution = fsolve(@(x) dH(x) - dH_st(i), x_initial_guess);
    x_st(i) = x_solution;
end

% Plotting Partial f_I Curve Mapping relationship 
figure(7)
scatter(x_st,dH_st);
hold on
xLine = linspace(min(x_st), max(x_st), 100);
yLine = dH(xLine);
plot(xLine, yLine, 'r-');
hold off
set(gca,'box','off');
xlabel('Estimated x_i^*');
ylabel('Estimated dH/dx_i^*');
xticks([0 0.55]);
yticks([0 270]);
alpha(0.8);
saveas(gcf,'PFIMapping.png');

% Plotting Estimated Input Current x* performance
figure(8)
scatter(x_star,x_st);
set(gca,'box','off');
xlabel('Ground Truth x_i^*');
ylabel('Estimated x_i^*');
xticks([0.3 0.55]);
yticks([0.3 0.55]);
alpha(0.8);
saveas(gcf,'x_starPerformance.png');

% Plotting f-I Curve Mapping relationship 
figure(9)
scatter(x_st,H(x_st));
hold on
xLine = linspace(min(x_st), max(x_st), 100);
yLine = H(xLine);
plot(xLine, yLine, 'r-');
hold off
set(gca,'box','off');
xlabel('Estimated x_i^*');
ylabel('Firing rate H');
xticks([0 0.55]);
yticks([0 3 10]);
alpha(0.8);
saveas(gcf,'FIMapping.png');

%% Evaluation of Outer Input I Estimation

I_recon = x_st - w_recon.*J.*S_star - G.*J.*C_recon*S_star;

figure(10)
scatter(I,I_recon);
hold on
plot([0.3,max(I)],[0.3,max(I)],'--','Color','k');
hold off
set(gca,'box','off');
xlabel('Ground Truth I_i');
ylabel('Estimated I_i');
xticks([0.3 0.35]);
yticks([0.3 0.35]);
alpha(0.8);
saveas(gcf,'IPerformance.png');

%% Testing the error contribution of components from estimated w, I

% Plotting Relative Error of w, I with Stable points S*
figure(11)
scatter(S_star,abs(w-w_recon)./w);
set(gca,'box','off');
xlabel('Stable Gating Variable S^*');
ylabel('Relative Error of w');
xticks([0 0.8]);
yticks([0 0.1 0.2]);
saveas(gcf,'ErrorW.png');

figure(12)
scatter(S_star,abs(I-I_recon)./I);
set(gca,'box','off');
xlabel('Stable Gating Variable S^*');
ylabel('Relative Error of I');
xticks([0 0.8]);
yticks([0 0.05 0.15]);
saveas(gcf,'ErrorI.png');

close all


%% heatmap
figure(13)
h1 = heatmap(SC,'Colormap',hot,'GridVisible','off');

figure(14)
h2 = heatmap(C_recon,'Colormap',hot,'GridVisible','off');
