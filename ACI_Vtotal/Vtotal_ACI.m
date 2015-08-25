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
        
%         nLargeKv = sum(kv>0.75);
%         meanKv = mean( kv(kv<=0.75) );
%         stdKv = std( kv(kv<=0.75) );
%         covKv = stdKv / meanKv;
        
        covKv =  mean(kv) / std(kv);
        kv( kv>0.75 ) = normrnd(0.75, 0.75*covKv, sum(kv>0.75), 1);
        ep_fe = kv .* ep_fu;
        
%         nLargeEp = sum( ep_fe>0.004 );
%         meanEp = mean( ep_fe(ep_fe<=0.004) );
%         stdEp = std( ep_fe(ep_fe<=0.004) );
%         covEp = 0.20;
        covEp = std(ep_fu) / mean(ep_fu);
        ep_fe( ep_fe>0.004 ) = normrnd(0.004, 0.004*covEp, nLargeEp, 1);
              
    case {'W'}
        ep_cov = std(ep_fu) / mean(ep_fu);
        ep_fe = normrnd( 0.004, 0.004*ep_cov, Nsim, 1);
        ep_fe( ep_fe<0 ) = 0;
        ep_fe(0.75*ep_fu<0.004) = 0.75*ep_fu(0.75*ep_fu<0.004);
        
%         ep_fe = 0.004 * ones(Nsim, 1);
%         nLarger = sum( 0.75*ep_fu>=ep_fe );
%         meanSmaller = mean( ep_fu(0.75*ep_fu<ep_fe) );
%         stdSmaller = std( ep_fu(0.75*ep_fu<ep_fe) );
%         covSmaller = stdSmaller / meanSmaller;  
%         index1 = 0.75*ep_fu<ep_fe;
%         index2 = 0.75*ep_fu>=ep_fe;
%         ep_fe( index1 ) = 0.75*ep_fu(index1);
%         ep_fe( index2 ) = normrnd( 0.004, 0.004*covSmaller, nLarger, 1);
end

f_fe = ep_fe .* E_frp;
V_frp = 2*t_frp.*w_frp.*f_fe.*(sin(beta)+cos(beta)).*dfv ./ s_frp /1000;
V_frp( V_frp<0 ) = 0;

% ACI 318-11: Steel Contribution

nLarger = sum( f_s>420 );
% meanSmaller = mean( f_s(f_s<=420) );
% stdSmaller = std( f_s(f_s<=420) );
% covSmaller = stdSmaller/meanSmaller;
covLarger = 0.10;
f_s( f_s>420 ) = normrnd(420, 420*covLarger, nLarger,1);

Av = 2*pi/4*D_bar.^2;
Vs = Av .* f_s .* d ./ s_bar /1000;
Vs( s_bar == 0 ) = 0;
Vs( Vs<0 ) = 0;

% contribution of concrete

Av_min = 0.062*sqrt(f_c) .* b .* s_bar ./ f_s;
temp = (0.35*b.*s_bar) ./ f_s;
Av_min( Av_min < temp ) = temp(Av_min < temp); 
sqrt_fc = sqrt(f_c);
isAdjust = (Av<Av_min) & (sqrt_fc>8.3);
nAdjust = sum( (Av<Av_min) & (sqrt_fc>8.3) );
covAdjust = 0.2;
sqrt_fc( isAdjust ) = normrnd( 8.3, 8.3*covAdjust, nAdjust, 1);

Vc = 0.17*sqrt_fc.*b.*d/1000;
Vc(Vc<0) = 0;

% over reinforcement check

Vr = Vs + V_frp;
temp = Vr>0.66*sqrt_fc.*b.*d/1000;
Vr( temp ) = 0.66*sqrt_fc(temp).*b(temp).*d(temp)/1000;

Vtotal = Vc + Vr;

return
end


