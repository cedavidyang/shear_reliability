function [Vtotal V_frp_interaction Vs Vc]= Vtotal_HK_design(db, FRP_tech, gamma_frpb, gamma_frp)
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

Nsim = 1;

% mean values to nomnal values
vc = 0.2;
f_cu = f_c ./ 0.8;
f_cu = f_cu * ( 1-1.645*vc );
vfrp=0.15;
f_frp = f_frp * (1-1.645*vfrp);
f_s = f_s ./ 1.13;

% partial safety factors
gamma_c = 1.25;
gamma_s = 1.15;

zt = dfrpt;
zb = (d-(h-dfrp))-0.1*d;
h_frp_e = zb-zt;
row = Nsim;

beta = beta/180*pi;
ratio = w_frp ./ (s_frp.*sin(beta));
temp = (2-ratio)./(1+ratio);
temp( temp<0 ) = 0;
beta_w = sqrt(temp);

if strcmp(FRP_tech, 'W')
    l_max = -1;
elseif strcmp(FRP_tech, 'U')
    l_max = h_frp_e ./ sin(beta);
elseif strcmp(FRP_tech, 'Side')
    l_max = h_frp_e ./ 2 ./ sin(beta);
end

% no partial safety factor is used for f_c in the determinatin of Vfrp
l_e = sqrt( E_frp.*t_frp ./ sqrt(0.8*f_cu) );
lambda = l_max ./ l_e;

beta_l = zeros(row, 1);
beta_l( lambda >= 1 ) = 1;
beta_l( lambda < 1 ) = sin(pi*lambda(lambda < 1)/2);

% Partial safety factor for debonding strength
s_debonding = 0.315*beta_w.*beta_l.*sqrt(E_frp.*sqrt(0.8*f_cu)./t_frp);
s_debonding = s_debonding / gamma_frpb;

% gamma_frp = 1.25 is used for FRP rupture failure
f_frp = f_frp/gamma_frp;

fail_mode = zeros(row, 1);
% 0.8f_frp is used to account for FRP strength degradation in shear
% strengthening
if strcmp(FRP_tech, 'W')
    s_frp_max = 0.8*f_frp;
    fail_mode = ones(Nsim,1); % 1 for fiber rupture
else
    s_frp_max = transpose( min( [s_debonding'; f_frp'*0.8] ) );
    fail_mode( s_debonding > 0.8*f_frp ) = 1;
    fail_mode( s_debonding <= 0.8*f_frp ) = 2; % 2 for debonding failure
end

D_frp = zeros(row, 1);
zeta = zt./zb;
D_frp( fail_mode==1 ) = (1+zeta(fail_mode==1))/2;
D_frp( (fail_mode==2 & lambda>1) ) = 1 - (pi-2)./(pi*lambda( (fail_mode==2 & lambda>1) ) );
D_frp( (fail_mode==2 & lambda<=1) ) = 2/pi./lambda((fail_mode==2 & lambda<=1)) .* (1-cos(pi/2*lambda((fail_mode==2 & lambda<=1)))) ./ sin(pi/2*lambda((fail_mode==2 & lambda<=1)));

f_frp_e = s_frp_max .* D_frp;
V_frp = 2*f_frp_e.*t_frp.*w_frp.*h_frp_e.*(sin(beta)+cos(beta))./s_frp/1000;
V_frp( V_frp < 0 ) = 0;

% steel-FRP interaction
if strcmp(FRP_tech, 'Side')
    % nominal value of steel yield strength is used in determining phi_s
    phi_s = zeros(Nsim, 1); A=zeros(Nsim, 1); Vs = zeros(Nsim, 1);
    mu = zeros(Nsim, 1); kfrp = zeros(Nsim, 1);
    
    phi_s( bar_type == 0 ) = 1e5 ./ (D_bar( bar_type == 0 ).^1.13 .* f_s( bar_type == 0 ).^ 1.71);
    A( bar_type == 0 ) = phi_s( bar_type == 0 ) .* (2.045*2*sin(beta( bar_type == 0 )).*lambda( bar_type == 0 )+3.24);
    phi_s( bar_type == 1 ) = 1e5 ./ (D_bar( bar_type == 1 ).^0.834 .* f_s( bar_type == 1 ).^ 1.88);
    A( bar_type == 1 ) = phi_s( bar_type == 1 ) .* (1.01*2*sin(beta( bar_type == 1 )).*lambda( bar_type == 1 )+2.13);
    
    Vs( bar_type ~= 2 ) = 1e-3*2*pi*D_bar( bar_type ~= 2 ).^2/4 .*f_s( bar_type ~= 2 ) .* d( bar_type ~= 2 ) ./ s_bar( bar_type ~= 2 );
    Vs( bar_type == 2 )  =0;
    mu( bar_type ~= 2 ) = Vs( bar_type ~= 2 ) ./ V_frp( bar_type ~= 2 ) ;
    mu( V_frp == 0 ) = 0;
    kfrp( bar_type ~= 2 ) = A( bar_type ~= 2 )./(A( bar_type ~= 2 )+mu( bar_type ~= 2 ));
    kfrp( bar_type == 2 ) = 1;
else
    f_s = f_s/gamma_s;
    Vs = 1e-3*2*pi*D_bar.^2/4 .*f_s .* d ./ s_bar;
    Vs( s_bar == 0 ) = 0;
    Vs( Vs<0 ) = 0;
    kfrp = 1;
end
    
V_frp_interaction = 2*kfrp.*f_frp_e.*t_frp.*w_frp.*h_frp_e.*(sin(beta)+cos(beta))./s_frp/1000;
V_frp_interaction( V_frp_interaction<0 ) = 0;

% concrete contribution
% partial safety factor for concrete: gamma_c = 1.25

temp1 = (400./ d).^(1/4);
vr = 0.4*ones(length(temp1),1);
vr( f_cu>40 & f_cu<80 ) = 0.4*(f_cu( f_cu>40 & f_cu<80 )/40).^(2/3);
vr( f_cu>=80 ) = 0.4*(80/40).^(2/3);
As_min = vr.* b .* s_bar ./ (0.87*f_s);
As = pi/4* D_bar.^2;
temp1( As >= As_min & temp1<1 ) = 1;

Vc = 0.79*(100*0.01).^(1/3).*temp1.*b.*d*1e-3;

% f_cu is nominal value
Vc( f_cu>25 & f_cu<80 ) = Vc( f_cu>25 & f_cu<80 ).*(f_cu( f_cu>25 & f_cu<80 )/25).^(1/3);
Vc( f_cu>=80 ) = Vc( f_cu>=80 ).*(80/25).^(1/3);
Vc = Vc / gamma_c;
Vc( Vc<0 ) = 0;

Vtotal = V_frp_interaction + Vs + Vc;
temp3 = min(0.8*sqrt(f_cu), 7);
Vtotal( Vtotal*1000./(b.*d)>temp3 ) = temp3(Vtotal*1000./(b.*d)>temp3).*b(Vtotal*1000./(b.*d)>temp3).*d(Vtotal*1000./(b.*d)>temp3)/1000;

return
end