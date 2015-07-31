% construction of design cases
clear; clc;
% beams
h = [200; 600];
as = [35; 60];
h2b = 2;
% L = [];

fc = [15; 70];

fs = [235; 400];
bar_type = [0; 1];
sd = [6; 12]; % reabr diameter
ss_min = 70;
ss_max = [150; 250];
% ss_max is related to h
s2d = [1; 3];
% 0 for deformed rebars
% 1 for round rebars

beta = [90; 45];
FRP_type = [0; 1];
% 0 for CFRP
% 1 for GFRP
Ef_C = [100e3; 300e3];
ef_C = [1/100; 1.5/100];
Ef_G = [20e3; 50e3];
ef_G = [2/100; 4/100];
tf = [0.2; 0.8];
str_type = [0; 1]; 
% 0 for continuous strengthening
% 1 for discrete strengthening

db_design = zeros(2^13, 19);
i = 1;

fc_db = zeros(2^13, 1);
b_db = zeros(2^13, 1);
h_db = zeros(2^13, 1);
d_db = zeros(2^13, 1);
dfrp_db = zeros(2^13, 1);
dfrpt_db = zeros(2^13, 1);
s2d_db = zeros(2^13, 1);
bar_type_db = zeros(2^13, 1);
sd_db = zeros(2^13, 1);
ss_db = zeros(2^13, 1);
fs_db = zeros(2^13, 1);
beta_db = zeros(2^13, 1);
FRP_type_db = zeros(2^13, 1);
str_type_db = zeros(2^13, 1);
Ef_db = zeros(2^13, 1);
tf_db = zeros(2^13, 1);
ff_db = zeros(2^13, 1);
wf_db = zeros(2^13, 1);
sf_db = zeros(2^13, 1);

for i_fc = 1:2
    for i_h = 1:2
        for i_s2d = 1:2
            for i_bar = 1:2
                for i_sd = 1:2
                    for i_ss = 1:2
                        for i_fs = 1:2
                            for i_beta = 1:2
                                for i_FRP = 1:2
                                    for i_Ef = 1:2
                                        for i_ef = 1:2
                                            for i_tf = 1:2
                                                for i_str = 1:2
                                                    fc_db(i) = fc(i_fc);
                                                    h_db(i) = h(i_h);
                                                    b_db(i) = h(i_h) / h2b;
                                                    d_db(i) = (h_db(i)==h(1))*(h_db(i)-as(1)) + (h_db(i)==h(2))*(h_db(i) - as(2));
                                                    dfrp_db(i) = h_db(i);
                                                    s2d_db(i) = s2d(i_s2d);
                                                    bar_type_db(i) = bar_type(i_bar);
                                                    sd_db(i) = sd(i_sd);
                                                    ss_db(i) = (i_ss==1)*ss_min + (i_ss==2 && h_db(i) == h(1)) * ss_max(1) + (i_ss==2 && h_db(i) == h(2)) * ss_max(2);
                                                    fs_db(i) = fs(i_fs);
                                                    beta_db(i) = beta(i_beta);
                                                    FRP_type_db(i) = FRP_type(i_FRP);
                                                    str_type_db(i) = str_type(i_str);
                                                    Ef_db(i) = (i_FRP==1)*Ef_C(i_Ef) + (i_FRP==2)*Ef_G(i_Ef);
                                                    tf_db(i) = tf(i_tf);
                                                    str_type_db(i) = str_type(i_str);
                                                    ff_db(i) = Ef_db(i) * ( (i_FRP==1)*ef_C(i_ef) + (i_FRP==2)*ef_G(i_ef) );
                                                    wf_db(i) = (i_str==1)*1 + (i_str==2)*50;
                                                    sf_db(i) = (i_str==1)*1/sin(beta(i_beta)/180*pi) + (i_str==2)*100;
                                                    i = i+1;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

db_design(:, 1) = fc_db;
db_design(:, 2) = b_db;
db_design(:, 3) = h_db;
db_design(:, 4) = d_db;
db_design(:, 5) = dfrp_db;
db_design(:, 6) = dfrpt_db;
db_design(:, 7) = s2d_db;
db_design(:, 8) = bar_type_db;
db_design(:, 9) = sd_db;
db_design(:, 10) = ss_db;
db_design(:, 11) = fs_db;
db_design(:, 12) = beta_db;
db_design(:, 13) = FRP_type_db;
db_design(:, 14) = str_type_db;
db_design(:, 15) = Ef_db;
db_design(:, 16) = tf_db;
db_design(:, 17) = ff_db;
db_design(:, 18) = wf_db;
db_design(:, 19) = sf_db;

db_design_HK = db_design;
save('db_design_HK.mat', 'db_design_HK')


                                                                    

