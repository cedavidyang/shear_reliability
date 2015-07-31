% size effects of HK model
clear; clc
load db_Side
load db_U
load db_W

[Vpre_Side, ~, ~, ~] = Vtotal_HK(db_Side, 'Side');
Vexp_Side = db_Side(:, 21);
ME_Side = Vexp_Side ./ Vpre_Side;

[Vpre_U, ~, ~, ~] = Vtotal_HK(db_U, 'U');
Vexp_U = db_U(:, 21);
ME_U = Vexp_U ./ Vpre_U;

[Vpre_W, ~, ~, ~] = Vtotal_HK(db_W, 'W');
Vexp_W = db_W(:, 21);
ME_W = Vexp_W ./ Vpre_W;

% db_U = []; ME_U = [];
db_Side = []; ME_Side = [];
db = [db_Side; db_U];
ME = [ME_Side; ME_U];

f_c  = db(:, 1);
b  = db(:, 2);
h  = db(:, 3);
d  = db(:, 4);
dfrp  = db(:, 5);
dfrpt  = db(:, 6);
s2d  = db(:, 7);
bar_type  = db(:, 8);
D_bar  = db(:, 9);
s_bar  = db(:, 10);
f_s  = db(:, 11);
beta  = db(:, 12);
FRP_type  = db(:, 13);
str_type  = db(:, 14);
E_frp  = db(:, 15);
t_frp  = db(:, 16);
f_frp  = db(:, 17);
w_frp  = db(:, 18);
s_frp  = db(:, 19);

% ME( D_bar~=0 ) = [];
% db( D_bar~=0, : ) = [];

ME( D_bar~= 0 | s2d<2.9 | s2d>3.1 ) = [];
db( D_bar~= 0 | s2d<2.9 | s2d>3.1, : ) = [];
d  = db(:, 4);

nbins = 9;
[n, ctr] = hist(d, nbins);
n = n';
d_lb = zeros(length(d), 1);

for ibins = 1:nbins
    if ibins == 1
        d_lb(1:n(1)) = ctr(ibins);
    else
        d_lb( sum(n(1:ibins-1))+1:sum(n(1:ibins)) ) = ctr(ibins);
    end
end

ME_mean = zeros(nbins, 1);
ME_std = zeros(nbins, 1);

for ibins = 1:nbins
    if ibins == 1
        ME_mean(1) = mean( ME(1:n(1)) );
        ME_std(1) = std( ME(1:n(1)) );     
    else
        ME_mean(ibins) = mean( ME( sum(n(1:ibins-1))+1:sum(n(1:ibins)) ) );
        ME_std(ibins) = std( ME( sum(n(1:ibins-1))+1:sum(n(1:ibins)) ) );
    end
end

% plot(d_lb, ME, 'o'); hold on;
% ctr(n==0) = []; ctr = ctr';
% ME_mean(n==0) = [];
% ME_std(n==0) = [];
% errorbar(ctr, ME_mean, ME_std, 'rs', 'markerface', 'r')
% xlabel('Effective depth (mm)'); ylabel('Model error');
% % xlim([50, 450]);
% xlim([200, 550]);
% 
% cf = fit(ctr, ME_mean, fittype('poly1'));
% yfit = cf.p1*ctr + cf.p2;
% y = ME_mean;
% R2 = norm(yfit -mean(y))^2/norm(y - mean(y))^2;
% hline = refline(cf.p1, cf.p2);
% set(hline, 'color', 'g', 'linewidth', 2);
% annotation(gcf, 'textbox', [0.59, 0.78, 0.29, 0.10], 'String', sprintf('y = %.5f x + %.5f\nR^2=%.5f', cf.p1, cf.p2, R2)); 

plot(d, ME, 'o'); hold on;
xlim([200, 600]);

cf = fit(d, ME, fittype('poly1'));
yfit = cf.p1*ctr + cf.p2;
y = ME;
R2 = norm(yfit -mean(y))^2/norm(y - mean(y))^2;
hline = refline(cf.p1, cf.p2);
set(hline, 'color', 'g', 'linewidth', 2);
annotation(gcf, 'textbox', [0.59, 0.78, 0.29, 0.10], 'String', sprintf('y= %.5f x + %.5f\nR^2=%.5f', cf.p1, cf.p2, R2)); 

xlabel('Effective depth (mm)')
ylabel('Model error')