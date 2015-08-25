% Reliability of U
clear; clc;
% Model Error
load db_U_ACI
[Vpre, ~, ~, ~] = Vtotal_ACI_ME(db_U, 'U');
Vexp = db_U(:, 21);
ModelError = Vexp ./ Vpre;
MEparam = lognfit(ModelError, 0.05);

% Monte Carlo Simulation of Resistance
IndexHat = 3.5; % target reliability index
factor = 0.5:0.05:1.00;
nFactor = length(factor);
factor = factor';
LD=transpose(0.25:0.25:3.0);
nLD = length(LD);

load db_design_ACI
db_design = db_design_ACI;
[nCase, ~] = size(db_design);
f_c  = db_design(:, 1);
b  = db_design(:, 2);
h  = db_design(:, 3);
d  = db_design(:, 4);
dfrp  = db_design(:, 5);
dfrpt  = db_design(:, 6);
s2d  = db_design(:, 7);
bar_type  = db_design(:, 8);
D_bar  = db_design(:, 9);
s_bar  = db_design(:, 10);
f_s  = db_design(:, 11);
beta  = db_design(:, 12);
FRP_type  = db_design(:, 13);
str_type  = db_design(:, 14);
E_frp  = db_design(:, 15);
t_frp  = db_design(:, 16);
f_frp  = db_design(:, 17);
w_frp  = db_design(:, 18);
s_frp  = db_design(:, 19);

nLC = nCase*nLD;
Nsim = 1e4;
RE = zeros(nLD, nFactor, nCase);
norm_RE = zeros(nFactor,1);
mean_RE = zeros(nFactor,1);
std_RE = zeros(nFactor,1);
upper_RE = zeros(nFactor,1);
lower_RE = zeros(nFactor,1);

matlabpool(6)
parfor i=1:nCase
    % nominal values
    b_nom = b(i);
    h_nom = h(i);
    d_nom = d(i);
    dfrp_nom = h_nom;
    dfrpt_nom = dfrpt(i);
    beta_nom = beta(i);
    t_frp_nom = t_frp(i);
    E_frp_nom = E_frp(i);
    f_frp_nom = f_frp(i);
    w_frp_nom = w_frp(i);
    s_frp_nom = s_frp(i);
    f_c_nom = f_c(i);
    D_bar_nom = D_bar(i);
    s_bar_nom = s_bar(i);
    f_s_nom = f_s(i);
    
    % sample values for MC simulation of resistance
    b_mean = b_nom + 2.54;
    b_std = 3.658;
    b_smp = normrnd( b_mean, b_std, Nsim, 1 );
    
    h_mean = h_nom - 3.05;
    h_std = 6.35;
    h_smp = normrnd(h_mean, h_std, Nsim, 1);
    
    d_mean = d_nom - 4.7;
    d_std = 12.7;
    d_smp = normrnd(d_mean, d_std, Nsim, 1);
    
    dfrp_smp = h_smp;
    dfrpt_smp = dfrpt_nom;
    
    beta_mean = beta_nom;
    beta_std = 1;
    beta_smp = normrnd(beta_mean, beta_std, Nsim, 1);
    
    t_frp_smp = t_frp_nom;
    
    E_frp_smp = E_frp_nom;
    
    f_frp_mean = f_frp_nom ./ (1-1.645*0.12);
    f_frp_std = 0.12*f_frp_mean;
    wblparam = fsolve(@(x) [x(1)*gamma(1+1./x(2)) - f_frp_mean ; x(1).^2 * (gamma(1+2./x(2)) - ( gamma(1+1./x(2)).^2)) - f_frp_std^2],[f_frp_mean;1.2/(f_frp_std/f_frp_mean)], optimset('Display','off'));
    f_frp_smp = wblrnd(wblparam(1), wblparam(2), Nsim, 1);
    
    w_frp_smp = w_frp_nom;
    s_frp_smp = s_frp_nom;
    
    f_c_mean = f_c_nom/(1-1.645*0.2);
    f_c_std = 0.2*f_c_mean;
%     f_c_mean = f_c_nom*1.25;
%     f_c_std = 0.2*f_c_mean;
    f_c_smp = normrnd(f_c_mean, f_c_std, Nsim, 1);
    
    D_bar_smp = D_bar_nom*ones(Nsim,1);
    s_bar_smp = s_bar_nom*ones(Nsim,1);
    
    f_s_mean = f_s_nom*1.13;
    f_s_std = 0.1*f_s_mean;
    f_s_log_std = sqrt( log( 0.1^2 + 1 ) );
    f_s_log_mean = log( f_s_mean ) - .5*f_s_log_std.^2;
    f_s_smp = lognrnd( f_s_log_mean, f_s_log_std, Nsim, 1);
        
    % database construction, sample values
    db_smp = zeros(Nsim, 19);
    db_smp(:, 1) = f_c_smp;
    db_smp(:, 2) = b_smp;
    db_smp(:, 3) = h_smp;
    db_smp(:, 4) = d_smp;
    db_smp(:, 5) = dfrp_smp;
    db_smp(:, 6) = dfrpt_smp;
    db_smp(:, 7) = s2d(i);
    db_smp(:, 8) = bar_type(i);
    db_smp(:, 9) = D_bar_smp;
    db_smp(:, 10) = s_bar_smp;
    db_smp(:, 11) = f_s_smp;
    db_smp(:, 12) = beta_smp;
    db_smp(:, 13) = FRP_type(i);
    db_smp(:, 14) = str_type(i);
    db_smp(:, 15) = E_frp_smp;
    db_smp(:, 16) = t_frp_smp;
    db_smp(:, 17) = f_frp_smp;
    db_smp(:, 18) = w_frp_smp;
    db_smp(:, 19) = s_frp_smp;
    db_smp( f_c_smp <= 0, : ) = [];
    
    % mean and nomonal values
    db_mean = zeros(1, 19);
    db_mean(:, 1) = f_c_mean;
    db_mean(:, 2) = b_nom;
    db_mean(:, 3) = h_nom;
    db_mean(:, 4) = d_nom;
    db_mean(:, 5) = dfrp_nom;
    db_mean(:, 6) = dfrpt_nom;
    db_mean(:, 7) = s2d(i);
    db_mean(:, 8) = bar_type(i);
    db_mean(:, 9) = D_bar_nom;
    db_mean(:, 10) = s_bar_nom;
    db_mean(:, 11) = f_s_mean;
    db_mean(:, 12) = beta_nom;
    db_mean(:, 13) = FRP_type(i);
    db_mean(:, 14) = str_type(i);
    db_mean(:, 15) = E_frp_nom;
    db_mean(:, 16) = t_frp_nom;
    db_mean(:, 17) = f_frp_mean;
    db_mean(:, 18) = w_frp_nom;
    db_mean(:, 19) = s_frp_nom;  
    
    [Vtotal, ~, ~, ~] = Vtotal_ACI(db_smp, 'U');

    n_warning = length( find(Vtotal<= 0) );
    if n_warning > 0
        fprintf('Vtotal <= 0, %d times\n', n_warning );
    end    
    Vtotal_design = zeros(nFactor, 1);
    RE_tmp = zeros(nLD, nFactor);
    for i_factor = 1:nFactor
        gamma_frp = 1.40;
        [Vtotal_design(i_factor),~,~,~] = Vtotal_ACI_design(db_mean, 'U', factor(i_factor));
        for i_LD=1:nLD
            Rparam = [mean(Vtotal), std(Vtotal)/mean(Vtotal), Vtotal_design(i_factor)];
%             Rparam = [mean(Vs+Vc), std(Vs+Vc)/mean(Vs+Vc), Vs_design+Vc_design];
            RE_tmp(i_LD, i_factor) = form_ACI(MEparam, LD(i_LD), Rparam );
        end
        RE(:, :, i) = RE_tmp;
    end
end
matlabpool close;

REdata_ACI_U = RE;
% save('REdata_ACI_U_12Frp_detConc.mat', 'REdata_ACI_U');
save('REdata_ACI_U_12Frp_detConc.mat', 'REdata_ACI_U');
RE_col = zeros(nLD, nCase);
for i_factor = 1:nFactor
    RE_col = RE(:, i_factor, :);
    RE_col = RE_col(:);
    RE_col( isnan(RE_col) ) = [];
    norm_RE(i_factor) = mean((RE_col-IndexHat).^2);
    mean_RE(i_factor) = mean(RE_col);
    std_RE(i_factor) = std(RE_col);
    upper_RE(i_factor) = prctile(RE_col, 95);
    lower_RE(i_factor) = prctile(RE_col, 5);
end

plot(factor, mean_RE, 'MarkerFaceColor','b','Marker','o',...
    'DisplayName','New model (Lognormal ME)', 'LineWidth', 2, 'Color', 'b', 'LineStyle', '-');
hold on
plot(factor, upper_RE, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
plot(factor, lower_RE, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
plot(factor, mean_RE+std_RE, 'LineWidth', 1, 'Color', 'g', 'LineStyle', '-.');
plot(factor, mean_RE-std_RE, 'LineWidth', 1, 'Color', 'g', 'LineStyle', '-.');
xlabel('Partial Safety Factor for FRP contribution');
ylabel('Reliability Index   \beta');
figure;
plot(factor, norm_RE, 'bo-', 'linewidth', 2, 'markerfacecolor', 'b')
xlabel('Partial Safety Factor for FRP contribution');
ylabel('Norm of Reliability Index');
