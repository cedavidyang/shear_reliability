function [Vtotal, Vfrp, Vs, Vc]= Vtotal_GB(db, FRP_tech)
% Determine the contribution of FRP, New HK guidelines
% Considering the steel-FRP interaction

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

f_cu = f_c / 0.8;

% GB notation
vc = 0.2;
f_cum = f_cu;
% f_t = 0.395 .* f_cum.^0.55 .* 0.8;
f_t = 0.395 .* f_cum.^0.55;

beta = beta / 180 *pi;
zt = dfrpt;
zb = (d-(h-dfrp))-0.1*d;
h_frp_e = zb-zt;
% [Nsim, ~] = size(db);

lambda = s2d;
lambda( lambda<1.5 ) = 1.5; lambda( lambda>3 ) = 3;
Vc = 1.75 ./ (lambda + 1) .* f_t .* b.*d;
As = pi*D_bar.^2/4;
Vs = f_s .* 2.*As./s_bar .* d;
Vs( s_bar==0 ) = 0;
Vcs = Vc + Vs;

if strcmp(FRP_tech, 'W')
%     phi0 = ones(Nsim, 1);
%     phi0( Vcs>0.7*f_t.*b.*d ) = 1 - ( Vcs - 0.7*f_t.*b.*d ) ./ ();
%    assume phi0 = 1
    phi0 = 1;
    lambda_Ef = 2*w_frp.*t_frp ./ ( b.*s_frp.*sin(beta) ) .* E_frp ./ f_t;
    e_fe_v = 8 ./ ( sqrt(lambda_Ef)+10 ) .* f_frp ./ E_frp;
    sf_vd = transpose( min([f_frp'; (E_frp.*e_fe_v)']) );
    Vfrp = phi0 * 2 .* w_frp .* t_frp ./ (s_frp.*sin(beta)) .* sf_vd .* h_frp_e .* ( sin(beta) + cos(beta) );
    Vtotal = Vfrp + Vcs;
else
    phi = (strcmp(FRP_tech, 'Side'))*1.0 + (strcmp(FRP_tech, 'U'))*1.3;
    beta_w = sqrt( (2.25-w_frp./(s_frp.*sin(beta)))./(1.25+w_frp./(s_frp.*sin(beta))) );
    tb = 1.2 * beta_w .* f_t;
    Kf = phi .* sin(beta) .* sqrt(E_frp.*t_frp) ./ ( sin(beta).*sqrt(E_frp.*t_frp) + 0.3*h_frp_e.*f_t ) ;
    Vfrp = Kf .* tb .* w_frp .* h_frp_e.^2 ./ (s_frp.*sin(beta)) .* ( sin(beta) + cos(beta) );
    Vtotal = Vfrp + Vcs;
end

% elimination of diagonal-compression failure
beta_c = 1.0;
Vtotal( Vtotal>0.2*beta_c.*f_c.*b.*d ) = 0.2 * beta_c .* f_c( Vtotal>0.2*beta_c.*f_c.*b.*d ) .* b( Vtotal>0.2*beta_c.*f_c.*b.*d ) .* d( Vtotal>0.2*beta_c.*f_c.*b.*d );
Vtotal = Vtotal*1e-3;
Vfrp = Vfrp*1e-3;
Vc = Vc*1e-3;
Vs = Vs*1e-3;

return
end