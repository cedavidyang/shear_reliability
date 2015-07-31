% extraction databases from xls files

%% Side bonding
[nu txt raw] = xlsread('.\Database\shear_database.xlsx', 'Side');
% geometric properties
b = nu(:,5); h=nu(:,6); d=nu(:,7); dfrp = nu(:,8); dfrpt = nu(:,9); s2d=nu(:, 12); 
row = length(h);
b = b(4:row); h = h(4:row); d=d(4:row); dfrp = dfrp(4:row); dfrpt = dfrpt(4:row);
s2d = s2d(4:row);
% strengthening scheme
FRP_type_txt = txt(4:row, 19);
n_txt = length(FRP_type_txt);
FRP_type = zeros(n_txt, 1);
for i_txt=1:n_txt
    FRP_type(i_txt) = ( strcmp(FRP_type_txt{i_txt},'CFRP') ) * 0 + ...
                      ( strcmp(FRP_type_txt{i_txt},'GFRP') ) * 1 + ...
                      ( ~strcmp(FRP_type_txt{i_txt},'CFRP') && ~strcmp(FRP_type_txt{i_txt},'GFRP')) * 2;
end
str_type = zeros(n_txt, 1);
% steel stirrups
bar_type = txt(4:row, 13);
D_bar = nu(4:row, 14);
s_bar = nu(4:row, 15);
f_s = nu(4:row, 17);
% other properties
beta = nu(4:row, 28); t_frp = nu(4:row, 23); E_frp = nu(4:row, 22);
f_frp = nu(4:row,24); w_frp = nu(4:row,25); s_frp = nu(4:row,26);
f_c = nu(4:row, 4); Vtotal_exp = nu(4:row, 30); Vf_exp = nu(4:row, 31);

D_bar( strcmp(bar_type, 'NA') ) = 0;
s_bar( strcmp(bar_type, 'NA') ) = 0;
f_s( strcmp(bar_type, 'NA') ) = 0;

bar_type_txt = bar_type;  bar_type = zeros(n_txt, 1);
for i_txt=1:n_txt
    bar_type(i_txt) = ( strcmp(bar_type_txt{i_txt},'D') ) * 0 + ...
                      ( strcmp(bar_type_txt{i_txt},'R') ) * 1 + ...
                      ( ~strcmp(bar_type_txt{i_txt},'D') && ~strcmp(bar_type_txt{i_txt},'R')) * 2;
end

row = length(dfrp);
db_Side = zeros(row, 21);
db_Side(:, 1) = f_c;
db_Side(:, 2) = b;
db_Side(:, 3) = h;
db_Side(:, 4) = d;
db_Side(:, 5) = dfrp;
db_Side(:, 6) = dfrpt;
db_Side(:, 7) = s2d;
db_Side(:, 8) = bar_type;
db_Side(:, 9) = D_bar;
db_Side(:, 10) = s_bar;
db_Side(:, 11) = f_s;
db_Side(:, 12) = beta;
db_Side(:, 13) = FRP_type;
db_Side(:, 14) = str_type;
db_Side(:, 15) = E_frp;
db_Side(:, 16) = t_frp;
db_Side(:, 17) = f_frp;
db_Side(:, 18) = w_frp;
db_Side(:, 19) = s_frp;
db_Side(:, 20) = Vf_exp;
db_Side(:, 21) = Vtotal_exp;

save('db_Side', 'db_Side');

%% U-jacketing
[nu txt raw] = xlsread('.\Database\shear_database.xlsx', 'U');
% geometric properties
b = nu(:,5); h=nu(:,6); d=nu(:,7); dfrp = nu(:,8); dfrpt = nu(:,9); s2d=nu(:, 12); 
row = length(h);
b = b(4:row); h = h(4:row); d=d(4:row); dfrp = dfrp(4:row); dfrpt = dfrpt(4:row);
s2d = s2d(4:row);
% strengthening scheme
FRP_type_txt = txt(4:row, 19);
n_txt = length(FRP_type_txt);
FRP_type = zeros(n_txt, 1);
for i_txt=1:n_txt
    FRP_type(i_txt) = ( strcmp(FRP_type_txt{i_txt},'CFRP') ) * 0 + ...
                      ( strcmp(FRP_type_txt{i_txt},'GFRP') ) * 1 + ...
                      ( ~strcmp(FRP_type_txt{i_txt},'CFRP') && ~strcmp(FRP_type_txt{i_txt},'GFRP')) * 2;
end
str_type = zeros(n_txt, 1);
% steel stirrups
bar_type = txt(4:row, 13);
D_bar = nu(4:row, 14);
s_bar = nu(4:row, 15);
f_s = nu(4:row, 17);
% other properties
beta = nu(4:row, 28); t_frp = nu(4:row, 23); E_frp = nu(4:row, 22);
f_frp = nu(4:row,24); w_frp = nu(4:row,25); s_frp = nu(4:row,26);
f_c = nu(4:row, 4); Vtotal_exp = nu(4:row, 30); Vf_exp = nu(4:row, 31);

D_bar( strcmp(bar_type, 'NA') ) = 0;
s_bar( strcmp(bar_type, 'NA') ) = 0;
f_s( strcmp(bar_type, 'NA') ) = 0;

bar_type_txt = bar_type;  bar_type = zeros(n_txt, 1);
for i_txt=1:n_txt
    FRP_type(i_txt) = ( strcmp(bar_type_txt{i_txt},'D') ) * 0 + ...
                      ( strcmp(bar_type_txt{i_txt},'R') ) * 1 + ...
                      ( ~strcmp(bar_type_txt{i_txt},'D') && ~strcmp(bar_type_txt{i_txt},'R')) * 2;
end

row = length(dfrp);
db_U = zeros(row, 21);
db_U(:, 1) = f_c;
db_U(:, 2) = b;
db_U(:, 3) = h;
db_U(:, 4) = d;
db_U(:, 5) = dfrp;
db_U(:, 6) = dfrpt;
db_U(:, 7) = s2d;
db_U(:, 8) = bar_type;
db_U(:, 9) = D_bar;
db_U(:, 10) = s_bar;
db_U(:, 11) = f_s;
db_U(:, 12) = beta;
db_U(:, 13) = FRP_type;
db_U(:, 14) = str_type;
db_U(:, 15) = E_frp;
db_U(:, 16) = t_frp;
db_U(:, 17) = f_frp;
db_U(:, 18) = w_frp;
db_U(:, 19) = s_frp;
db_U(:, 20) = Vf_exp;
db_U(:, 21) = Vtotal_exp;

save('db_U', 'db_U');

%% Wrapping scheme

[nu txt raw] = xlsread('.\Database\shear_database.xlsx', 'W');
% geometric properties
b = nu(:,5); h=nu(:,6); d=nu(:,7); dfrp = nu(:,8); dfrpt = nu(:,9); s2d=nu(:, 12); 
row = length(h);
b = b(4:row); h = h(4:row); d=d(4:row); dfrp = dfrp(4:row); dfrpt = dfrpt(4:row);
s2d = s2d(4:row);
% strengthening scheme
FRP_type_txt = txt(4:row, 19);
n_txt = length(FRP_type_txt);
FRP_type = zeros(n_txt, 1);
for i_txt=1:n_txt
    FRP_type(i_txt) = ( strcmp(FRP_type_txt{i_txt},'CFRP') ) * 0 + ...
                      ( strcmp(FRP_type_txt{i_txt},'GFRP') ) * 1 + ...
                      ( ~strcmp(FRP_type_txt{i_txt},'CFRP') && ~strcmp(FRP_type_txt{i_txt},'GFRP')) * 2;
end
str_type = zeros(n_txt, 1);
% steel stirrups
bar_type = txt(4:row, 13);
D_bar = nu(4:row, 14);
s_bar = nu(4:row, 15);
f_s = nu(4:row, 17);
% other properties
beta = nu(4:row, 28); t_frp = nu(4:row, 23); E_frp = nu(4:row, 22);
f_frp = nu(4:row,24); w_frp = nu(4:row,25); s_frp = nu(4:row,26);
f_c = nu(4:row, 4); Vtotal_exp = nu(4:row, 30); Vf_exp = nu(4:row, 31);

D_bar( strcmp(bar_type, 'NA') ) = 0;
s_bar( strcmp(bar_type, 'NA') ) = 0;
f_s( strcmp(bar_type, 'NA') ) = 0;

bar_type_txt = bar_type;  bar_type = zeros(n_txt, 1);
for i_txt=1:n_txt
    FRP_type(i_txt) = ( strcmp(bar_type_txt{i_txt},'D') ) * 0 + ...
                      ( strcmp(bar_type_txt{i_txt},'R') ) * 1 + ...
                      ( ~strcmp(bar_type_txt{i_txt},'D') && ~strcmp(bar_type_txt{i_txt},'R')) * 2;
end

row = length(dfrp);
db_W = zeros(row, 21);
db_W(:, 1) = f_c;
db_W(:, 2) = b;
db_W(:, 3) = h;
db_W(:, 4) = d;
db_W(:, 5) = dfrp;
db_W(:, 6) = dfrpt;
db_W(:, 7) = s2d;
db_W(:, 8) = bar_type;
db_W(:, 9) = D_bar;
db_W(:, 10) = s_bar;
db_W(:, 11) = f_s;
db_W(:, 12) = beta;
db_W(:, 13) = FRP_type;
db_W(:, 14) = str_type;
db_W(:, 15) = E_frp;
db_W(:, 16) = t_frp;
db_W(:, 17) = f_frp;
db_W(:, 18) = w_frp;
db_W(:, 19) = s_frp;
db_W(:, 20) = Vf_exp;
db_W(:, 21) = Vtotal_exp;

save('db_W', 'db_W');
