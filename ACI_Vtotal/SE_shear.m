% size effects of ACI model
clear; clc
load db_Side_ACI
load db_U_ACI
load db_W_ACI

[Vpre_Side, ~, ~, ~] = Vtotal_ACI(db_Side, 'Side');
Vexp_Side = db_Side(:, 21);
ME_Side = Vexp_Side ./ Vpre_Side;

[Vpre_U, ~, ~, ~] = Vtotal_ACI(db_U, 'U');
Vexp_U = db_U(:, 21);
ME_U = Vexp_U ./ Vpre_U;

[Vpre_W, ~, ~, ~] = Vtotal_ACI(db_W, 'W');
Vexp_W = db_W(:, 21);
ME_W = Vexp_W ./ Vpre_W;

% db_U = []; Vexp_U = [];
db_Side = []; Vexp_Side = [];
db = [db_Side; db_U];
Vtotal = [Vexp_Side; Vexp_U]*1e3;

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

v = Vtotal ./ (b.*d.*sqrt(f_c));

db( D_bar~= 0 | s2d<2.9 | s2d>3.1, : ) = [];
v( D_bar~= 0 | s2d<2.9 | s2d>3.1, : ) = [];
d( D_bar~= 0 | s2d<2.9 | s2d>3.1, : ) = [];

plot(d, v, 'o'); hold on;
xlim([200, 600]);

cf = fit(d, v, fittype('poly1'));
yfit = cf.p1*d + cf.p2;
y = v;
R2 = norm(yfit -mean(y))^2/norm(y - mean(y))^2;
hline = refline(cf.p1, cf.p2);
set(hline, 'color', 'g', 'linewidth', 2);
annotation(gcf, 'textbox', [0.59, 0.78, 0.29, 0.10], 'String', sprintf('y= %.5f x + %.5f\nR^2=%.5f', cf.p1, cf.p2, R2)); 

xlabel('Effective depth (mm)', 'interpreter', 'latex')
ylabel('$$V / bd \sqrt{f_c}$$', 'interpreter', 'latex')