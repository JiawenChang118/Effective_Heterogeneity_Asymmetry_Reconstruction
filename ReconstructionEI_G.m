%% Parameters setting

load('SC.mat');
SC = fln;
SC = SC./max(SC,[],'all');
N = length(SC);
T = 5000; %s;
dt = 0.001; %1ms;
t = linspace(0,T,T/dt);

tau = [0.1 , 0.01]; % E and I
gamma = 0.641;
sigma = 0.01; % nA, the noise amplitud.
w_E = 1; % Scale I_0
w_I = 0.7; % Scale I_0
J = 0.15; %J_NMDA
I_0 = 0.382;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Heterogeneity
w_EE = 0.21*linspace(0.7,1.3,N)';

w_EI = 0.15*linspace(0.7,1.3,N)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_IE = 2*ones(N,1); %nA, J_i.         % Can be adjusted to set firing rate ~ 3Hz.

w = struct();
w.w_E = w_E;w.w_I = w_I;w.w_IE = w_IE;w.w_EI = w_EI;w.w_EE = w_EE;

% FIC
% Setting f-I function H and dH

a_E = 310; % nC
b_E = 125; % Hz
d_E = 0.16; % s

a_I = 615;
b_I = 177;
d_I = 0.087;

H_E = @(x)dMFM_H(x,a_E,b_E,d_E); % f-I Curve for Excitatory Population
H_I = @(x)dMFM_H(x,a_I,b_I,d_I); % f-I Curve for Inhibitory Population

dH_E = @(x) - 310./(exp(20 - (248.*x)./5) - 1) - ...
    (248.*exp(20 - (248.*x)./5).*(310.*x - 125))./(5*(exp(20 - (248.*x)./5) - 1).^2);
dH_I = @(x) - 615./(exp(15399/1000 - (10701.*x)./200) - 1) - ...
    (10701.*exp(15399/1000 - (10701.*x)./200).*(615.*x - 177))./(200*(exp(15399/1000 - (10701.*x)./200) - 1).^2);

k = 20;
% k = 30; 
% repeat_time = 10; 
repeat_time = 1;
g_vector = linspace(0,2,2*k+1);
g_vector = g_vector(2:end);

% Setting Iterations
SSE_J_G = zeros(repeat_time,length(g_vector));
SSE_SC_G = zeros(repeat_time,length(g_vector));
Corr_J_G = zeros(repeat_time,length(g_vector));
Corr_SC_G = zeros(repeat_time,length(g_vector));
RE_h_G = zeros(repeat_time,length(g_vector));
Corr_re = zeros(repeat_time,length(g_vector));
[S_E_st_low,S_E_st_high,H_E_st_low,H_E_st_high] = deal(zeros(N,length(g_vector)));
Eigenvalues_G = zeros(2*N,repeat_time,length(g_vector));
l = 0;

%% Reduced Wong-Wang Model Simulation
tic
options = optimoptions('fsolve', 'Display', 'off');

for q = 1:length(g_vector)
    for repeat = 1:repeat_time
        l = l+1;
        G = g_vector(q);

        % FIC
        [w_IE] = FIC(G,w_EE,w_EI);
        w.w_IE = w_IE;
        
        [S_E, ~, S_I, ~] = EI_dMFM(SC, dt, T, w, G, sigma, H_E, H_I, tau);
        S_E = S_E(:,1e6:end);
        S_I = S_I(:,1e6:end);
        S_E_star = mean(S_E,2);
        S_I_star = mean(S_I,2);
        x_E_star = w_EE.*S_E_star + G.*J.*SC*S_E_star - w_IE.*S_I_star + w_E*I_0;
        x_I_star = w_EI.*S_E_star - S_I_star + w_I*I_0;
        H_E_star = H_E(x_E_star);
        H_I_star = H_I(x_I_star);
        
        [S_E_high, ~, S_I_high, ~] = EI_dMFM_high(SC, dt, T, w, G, sigma, H_E, H_I, tau);
        S_E_star_high = mean(S_E_high,2);
        S_I_star_high = mean(S_I_high,2);
        x_E_star_high = w_EE.*S_E_star_high + G.*J.*SC*S_E_star_high - w_IE.*S_I_star_high + w_E*I_0;
        x_I_star_high = w_EI.*S_E_star_high - S_I_star_high + w_I*I_0;
        H_E_star_high = H_E(x_E_star_high);
        H_I_star_high = H_I(x_I_star_high);
        
        S_E_st_low(:,q) = S_E_st_low(:,q)+S_E_star;
        S_E_st_high(:,q) = S_E_st_high(:,q)+S_E_star_high;
        H_E_st_low(:,q) = H_E_st_low(:,q) + H_E_star;
        H_E_st_high(:,q) = H_E_st_high(:,q) + H_E_star_high;
        
        % Jacobian Matrix Calculation
        
        Jacob_EE = zeros(N);
        Jacob_EI = zeros(N);
        Jacob_IE = zeros(N);
        Jacob_II = zeros(N);

        for i = 1:N
            for j = 1:N
                if i == j
                    Jacob_EE(i,j) = -1/tau(1)-gamma*H_E(x_E_star(i))+w_EE(i)*gamma*(1-S_E_star(i))*dH_E(x_E_star(i));
                    Jacob_EI(i,j) = w_EI(i)*dH_I(x_I_star(i));
                    Jacob_IE(i,j) = -w_IE(i)*(1-S_E_star(i))*gamma*dH_E(x_E_star(i));
                    Jacob_II(i,j) = -1/tau(2) - dH_I(x_I_star(i));
                else
                    Jacob_EE(i,j) = gamma*G*J*(1-S_E_star(i))*SC(i,j)*dH_E(x_E_star(i));
                end
            end
        end
        Jacob = [Jacob_EE Jacob_IE
            Jacob_EI Jacob_II];
        y_ana = 1./(gamma*G*J.*(1-S_E_star).*dH_E(x_E_star));
        Jacob_eff = Jacob_EE - Jacob_IE/Jacob_II*Jacob_EI;
        [V,D] = eig(Jacob,'vector');
        
        % Linear Reconstruct Jacobian Matrix with S
        Jacob_est = LinearReconst(S_E,dt); % Estimating Jacobian Matrix
        
        % Reconstruct spatial properties with Jacobian Matrix
        SC_sym = (SC + SC')/2;
        [y_st,C_recon] = RevealHHetero2(SC_sym,Jacob_est);
        C_recon = C_recon-diag(diag(C_recon));
        h_st = 1./y_st;
        
        J_resim = diag(diag(Jacob_est))+h_st.*C_recon;

        % Evaluation of FC Resimulations
        [S_re] = LinearResim(J_resim, dt, T, sigma);

        % Evaluation of Jacobian Matrix J and SC Estimation
        FC1 = corr(S_E'); FC2 = corr(S_re');
        Corr_resim = corr(reshape(FC1-eye(N),N^2,1),reshape(FC2-eye(N),N^2,1));


        Eigenvalues_G(:,repeat,q) = D;
        SSE_J = norm(Jacob_eff-Jacob_est)/norm(Jacob_eff); 
        A_vec = reshape(Jacob_eff-diag(diag(Jacob_eff)),N^2,1);
        B_vec = reshape(Jacob_est-diag(diag(Jacob_est)),N^2,1); 
        Corr_J = corr(A_vec,B_vec);

        SSE_SC = norm(SC-C_recon)/norm(SC); 
        A_vec = reshape(SC-diag(diag(SC)),N^2,1);
        B_vec = reshape(C_recon-diag(diag(C_recon)),N^2,1); 
        Corr_SC = corr(A_vec,B_vec);

        Corr_re(repeat,q) = Corr_resim;
        SSE_J_G(repeat,q) = SSE_J;
        SSE_SC_G(repeat,q) = SSE_SC;
        Corr_J_G(repeat,q) = Corr_J;
        Corr_SC_G(repeat,q) = Corr_SC;
        RE_h_G(repeat,q) = norm(y_ana-y_st)/norm(y_ana);
        
    end
    S_E_st_low(:,q) = S_E_st_low(:,q)/repeat;
    S_E_st_high(:,q) = S_E_st_high(:,q)/repeat;
    H_E_st_low(:,q) = H_E_st_low(:,q)/repeat;
    H_E_st_high(:,q) = H_E_st_high(:,q)/repeat;
end
toc

%% Results Plotting
Color1 = [33,49,80]./256;
Color2 = [199,35,54]./256;
g_vector = linspace(0,2,2*k+1);
g_vector = g_vector(2:end);

% figure(1)
% % plot(g_vector,SSE_J_G(1,:),'-o','Color',[199/256,210/256,219/256],'LineWidth',1,'MarkerEdgeColor',...
% %     [199/256,210/256,219/256],'MarkerFaceColor',[199/256,210/256,219/256]);
% % hold on
% % plot(g_vector,SSE_SC_G(1,:),'-o','Color',[243/256,206/256,196/256],'LineWidth',1,'MarkerEdgeColor',...
% %     [243/256,206/256,196/256],'MarkerFaceColor',[243/256,206/256,196/256]);
% % hold off
% plot(g_vector,mean(SSE_J_G,1),'Color',Color1,'LineWidth',2);
% hold on
% plot(g_vector,mean(SSE_SC_G,1),'Color',Color2,'LineWidth',2);
% patch([g_vector fliplr(g_vector)],[mean(SSE_SC_G,1)-std(SSE_SC_G,0,1) fliplr(mean(SSE_SC_G,1)+std(SSE_SC_G,0,1))],...
%     Color2,'edgealpha', '0', 'facealpha', '.2')
% patch([g_vector fliplr(g_vector)],[mean(SSE_J_G,1)-std(SSE_J_G,0,1) fliplr(mean(SSE_J_G,1)+std(SSE_J_G,0,1))],...
%     Color1,'edgealpha', '0', 'facealpha', '.2')
% hold off
% legend('J','SC');
% set(gca,'box','off');
% xticks([1 2]);
% yticks([0 0.5 1]);
% xlabel('G');
% % alpha(0.8);
% ylabel('Relative Error');

% figure(2)
% % plot(g_vector,Corr_J_G,'-o','MarkerEdgeColor',[155/256,190/256,219/256],...
% %             'MarkerFaceColor',[155/256,190/256,219/256]);
% % hold on
% % plot(g_vector,Corr_SC_G,'-o','MarkerEdgeColor',[243/256,169/256,147/256],...
% %             'MarkerFaceColor',[243/256,169/256,147/256])
% % hold off
% plot(g_vector,mean(Corr_J_G,1),'Color',Color1,'LineWidth',2);
% hold on
% plot(g_vector,mean(Corr_SC_G,1),'Color',Color2,'LineWidth',2);
% patch([g_vector fliplr(g_vector)],[mean(Corr_SC_G,1)-std(Corr_SC_G,0,1) fliplr(mean(Corr_SC_G,1)+std(Corr_SC_G,0,1))],...
%     Color2,'edgealpha', '0', 'facealpha', '.2')
% patch([g_vector fliplr(g_vector)],[mean(Corr_J_G,1)-std(Corr_J_G,0,1) fliplr(mean(Corr_J_G,1)+std(Corr_J_G,0,1))],...
%     Color1,'edgealpha', '0', 'facealpha', '.2')
% set(gca,'box','off');
% legend('J','SC');
% % alpha(0.8);
% xticks([1 2]);
% yticks([0 0.5 1]);
% ylim([0 1]);
% xlabel('G');
% ylabel('Corr');

figure;
plot(g_vector,mean(SSE_J_G,1),'Color',Color2,'LineWidth',3);
hold on
patch([g_vector fliplr(g_vector)],[mean(SSE_J_G,1)-std(SSE_J_G,0,1) fliplr(mean(SSE_J_G,1)+std(SSE_J_G,0,1))],...
    Color2,'edgealpha', '0', 'facealpha', '.2')
set(gca,'box','off', 'LineWidth', 3, 'FontName', 'Arial', 'FontSize', 20);
legend('J', 'FontName', 'Arial', 'FontSize', 20);
xticks([1 2]);
yticks([0 0.5 1]);
xlabel('G');
ylabel('Relative Error');


figure;
plot(g_vector,mean(SSE_SC_G,1),'Color',Color1,'LineWidth',3);
hold on
plot(g_vector,mean(RE_h_G,1),'Color',Color2,'LineWidth',3);
patch([g_vector fliplr(g_vector)],[mean(SSE_SC_G,1)-std(SSE_SC_G,0,1) fliplr(mean(SSE_SC_G,1)+std(SSE_SC_G,0,1))],...
    Color1,'edgealpha', '0', 'facealpha', '.2')
patch([g_vector fliplr(g_vector)],[mean(RE_h_G,1)-std(RE_h_G,0,1) fliplr(mean(RE_h_G,1)+std(RE_h_G,0,1))],...
    Color2,'edgealpha', '0', 'facealpha', '.2')
% xline(g_vector(18),'--','LineWidth',3)
% xline(g_vector(29),'--','LineWidth',3)
set(gca,'box','off', 'LineWidth', 3, 'FontName', 'Arial', 'FontSize', 20);
legend('SC','h', 'FontName', 'Arial', 'FontSize', 20);
xticks([1 2]);
yticks([0 0.5 1]);
xlabel('G');
ylabel('Relative Error');

figure;
plot(g_vector,mean(Corr_re,1),'Color',Color2,'LineWidth',3);
hold on
patch([g_vector fliplr(g_vector)],[mean(Corr_re,1)-std(Corr_re,0,1) fliplr(mean(Corr_re,1)+std(Corr_re,0,1))],...
    Color2,'edgealpha', '0', 'facealpha', '.2');
% xline(0.5,'--','LineWidth',3)
% xline(g_vector(26),'--','LineWidth',3)
% xline(g_vector(30),'--','LineWidth',3)
hold off
set(gca,'box','off', 'LineWidth', 3, 'FontName', 'Arial', 'FontSize', 20);
xticks([1 2]);
yticks([0 0.5 1]);
xlabel('G');
ylabel('Correlation of FCs');


figure;
g_vector = g_vector(:);
x_low = repelem(g_vector, size(S_E_st_low, 1));
y_low = S_E_st_low(:);
x_high = repelem(g_vector, size(S_E_st_high, 1));
y_high = S_E_st_high(:);
scatter(x_low, y_low, 'filled', ...
    'MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', Color2, ...
    'MarkerFaceAlpha', 0.7);
hold on;
scatter(x_high, y_high, 'filled', ...
    'MarkerEdgeColor', 'none', ... 
    'MarkerFaceColor', Color1, ...
    'MarkerFaceAlpha', 0.7);
% xline(0.65,'--','LineWidth',3)
% xline(g_vector(26),'--','LineWidth',3)
% xline(g_vector(30),'--','LineWidth',3)
hold off;
set(gca, 'box', 'off', 'LineWidth', 3, 'FontName', 'Arial', 'FontSize', 20');
xticks([1 2]);
yticks([0 0.5 1]);
xlabel('G');
ylabel('S^E_{st}');

figure;
g_vector = g_vector(:);
x_low = repelem(g_vector, size(H_E_st_low, 1));
y_low = H_E_st_low(:);
x_high = repelem(g_vector, size(H_E_st_high, 1));
y_high = H_E_st_high(:);
scatter(x_low, y_low, 'filled', ...
    'MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', Color2, ...
    'MarkerFaceAlpha', 0.7);
hold on;
scatter(x_high, y_high, 'filled', ...
    'MarkerEdgeColor', 'none', ... 
    'MarkerFaceColor', Color1, ...
    'MarkerFaceAlpha', 0.7);
% xline(0.65,'--','LineWidth',3)
% xline(g_vector(26),'--','LineWidth',3)
% xline(g_vector(30),'--','LineWidth',3)
hold off;
set(gca, 'box', 'off', 'LineWidth', 3, 'FontName', 'Arial', 'FontSize', 20');
xticks([1 2]);
% yticks([0 0.5 1]);
xlabel('G');
ylabel('Regional Firing Rate');

% figure;
% plot(g_vector,mean(RE_h_G,1),'Color',Color1,'LineWidth',2);
% hold on
% patch([g_vector fliplr(g_vector)],[mean(RE_h_G,1)-std(RE_h_G,0,1) fliplr(mean(RE_h_G,1)+std(RE_h_G,0,1))],...
%     Color1,'edgealpha', '0', 'facealpha', '.2')
% set(gca,'box','off');
% legend('h');
% % alpha(0.8);
% xticks([1 2]);
% yticks([0 0.1 0.5]);
% % ylim([0 1]);
% xlabel('G');
% ylabel('Relative Error');
% set(gca, 'FontName', 'Arial')
