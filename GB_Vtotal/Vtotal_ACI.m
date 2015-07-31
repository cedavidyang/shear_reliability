function [Vtotal, V_frp, Vs, Vc] = Vtotal_ACI(db, FRP_tech)
% Determine the contribution of FRP, ACI's model (2008)

% reading database

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
[Nsim, ~] = size(db);

% ACI 440.2R: FRP contribution

dfv = d-dfrpt;
ep_fu = f_frp ./ E_frp;
beta = beta/180*pi;

switch FRP_tech
    case {'U', 'Side'}
        Le = 23300 ./ (t_frp .* E_frp).^0.58;
        k1 = (f_c ./ 27) .^ (2/3);
        if strcmp(FRP_tech, 'U')
            k2 = (dfv-Le) ./ dfv;
        elseif strcmp(FRP_tech, 'Side')
            k2 = (dfv-2*Le) ./ dfv;
        end
        
        k2( k2<0 ) = 0;
        kv = k1.*k2.*Le ./ (11900*ep_fu);
        kv( kv>0.75 ) = 0.75;
        
        ep_fe = kv .* ep_fu;
        ep_fe( ep_fe>0.004 ) = 0.004;
              
    case {'W'}
        ep_fe = 0.004 * ones(Nsim, 1);
        ep_fe( ep_fe>0.75*ep_fu ) = 0.75*ep_fu(ep_fe>0.75*ep_fu);
end

f_fe = ep_fe .* E_frp;
V_frp = 2*t_frp.*w_frp.*f_fe.*(sin(beta)+cos(beta)).*dfv ./ s_frp /1000;
V_frp( V_frp<0 ) = 0;

% ACI 318-11: Steel Contribution

f_s( f_s>420 ) = 420;
Av = 2*pi/4*D_bar.^2;
Vs = Av .* f_s .* d ./ s_bar /1000;
Vs( s_bar == 0 ) = 0;
Vs( Vs<0 ) = 0;

% contribution of concrete

Av_min = 0.062*sqrt(f_c) .* b .* s_bar ./ f_s;
temp = (0.35*b.*s_bar) ./ f_s;
Av_min( Av_min < temp ) = temp(Av_min < temp); 
sqrt_fc = sqrt(f_c);
sqrt_fc( Av<Av_min & sqrt_fc>8.3 ) = 8.3;

Vc = 0.17*sqrt_fc.*b.*d/1000;
Vc(Vc<0) = 0;

% over reinforcement check

Vr = Vs + V_frp;
temp = Vr>0.66*sqrt_fc.*b.*d;
Vr( temp ) = 0.66*sqrt_fc(temp).*b(temp).*d(temp);

Vtotal = Vc + Vr;

return
end


