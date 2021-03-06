function index = form_ACI(MEparam, LD, Rparam)
% return reliabilty index according to partial safty coefficient(PSC)
% MEparam: parameters of model error
% LD: ratio of Live Load to Dead Load
% Rparam: Resistance parameters
% index: reliability index
 
cov_ME = sqrt( exp( MEparam(2).^2 )-1 );
mu_ME = exp( MEparam(1)+0.5*MEparam(2).^2 );
std_ME = mu_ME*cov_ME;
mu_R = Rparam(1);
std_R = Rparam(2)*mu_R;
design_R = Rparam(3);

S = design_R;
L_norm = LD*S/(1.4+1.7*LD);
% mu_L = 1.2*L_norm;
mu_L = 1.0*L_norm;
sigma_L = 0.25*mu_L;
D_norm = S/(1.4+1.7*LD);
mu_D = 1.05*D_norm;
sigma_D = 0.10*mu_D;

%% FORM independent variables
% muX = [mu_ME; mu_R; mu_D; mu_L];
% sigmaX = [std_ME; std_R; sigma_D; sigma_L];
% aEv = sqrt(6)*sigmaX(4)/pi; uEv = -psi(1)*aEv-muX(4);
% muX1 = muX; sigmaX1 = sigmaX;
% x=muX; normX = eps;
% while abs(norm(x)-normX)/normX > 1e-6
%     normX = norm(x);
%     g = x(1)*x(2)-(x(3)+x(4));
%     gX = [x(2); x(1); -1; -1];
%     % equivalent normal
%     cdfX = [logncdf(x(1), MEparam(1), MEparam(2));
%             1-evcdf(-x(4), uEv, aEv)];
%     pdfX = [lognpdf(x(1), MEparam(1), MEparam(2));
%             evpdf(-x(4), uEv, aEv)];
%     nc = norminv(cdfX);
%     sigmaX1([1,4]) = normpdf(nc) ./ pdfX;
%     muX1([1,4]) = x([1,4]) - nc.*sigmaX1([1,4]);
%     gs = gX.*sigmaX1; alphaX = -gs/norm(gs);
%     index = (g+gX'*(muX1-x))/norm(gs);
%     x = muX1+index*sigmaX1.*alphaX;
% end

%% No consideration of model error
% ME_ptile = 1;
% muX = [mu_R; mu_D; mu_L];
% sigmaX = [std_R; sigma_D; sigma_L];
% aEv = sqrt(6)*sigmaX(3)/pi; uEv = -psi(1)*aEv-muX(3);
% muX1 = muX; sigmaX1 = sigmaX;
% x=muX; normX = eps;
% while abs(norm(x)-normX)/normX > 1e-6
%     normX = norm(x);
%     g = ME_ptile*x(1)-(x(2)+x(3));
%     gX = [ME_ptile; -1; -1];
%     % equivalent normal
%     cdfX = 1-evcdf(-x(3), uEv, aEv);
% %     n_cdfX = length(cdfX);
% %     for i_cdfX=1:n_cdfX
% %         if abs(cdfX(i_cdfX) - 0) <= eps
% %             cdfX = 0.001;
% %         elseif abs(cdfX(i_cdfX) - 1) <= eps
% %             cdfX = 0.999;
% %         end
% %     end
%     pdfX = evpdf(-x(3), uEv, aEv);
%     nc = norminv(cdfX);
% %     n_pdfX = length(pdfX);
% %     for i_pdfX=1:n_pdfX
% %         if abs(pdfX(i_pdfX) - 0) <= eps
% %             sigmaX1(3) = sigmaX(3);
% %         else
%             sigmaX1(3) = normpdf(nc) ./ pdfX;
% %         end
% %     end
%     muX1(3) = x(3) - nc.*sigmaX1(3);
%     gs = gX.*sigmaX1; alphaX = -gs/norm(gs);
%     index = (g+gX'*(muX1-x))/norm(gs);
%     x = muX1+index*sigmaX1.*alphaX;
% end

%% D and R is correlated
rho_RD = 0;
% calculate Nataf coefficient R, since D and R are normal distribution, R=1
muX = [mu_ME; mu_R; mu_D; mu_L];
sigmaX = [std_ME; std_R; sigma_D; sigma_L];
aEv = sqrt(6)*sigmaX(4)/pi; uEv = -psi(1)*aEv-muX(4);
muX1 = muX; sigmaX1 = sigmaX; index = 10; deltaBeta = 10;
rhoX1 = [1 0 0 0; 0 1 rho_RD 0; 0 rho_RD 1 0; 0 0 0 1];
x = muX; normX = eps;
while abs(norm(x)-normX)/normX>1e-8 || deltaBeta > 2e-3
    normX = norm(x);
    g = x(1)*x(2)-(x(3)+x(4));
    gX = [x(2); x(1); -1; -1];
    cdfX = [logncdf(x(1), MEparam(1), MEparam(2));
            normcdf(x(2), mu_R, std_R);
            normcdf(x(3), mu_D, sigma_D);
            1-evcdf(-x(4), uEv, aEv)];
    pdfX = [lognpdf(x(1), MEparam(1), MEparam(2));
            normpdf(x(2), mu_R, std_R);
            normpdf(x(3), mu_D, sigma_D);        
            evpdf(-x(4), uEv, aEv)];  
    nc = norminv(cdfX);
    sigmaX1 = normpdf(nc)./pdfX;
    muX1 = x - nc .* sigmaX1;
    gs = gX.*sigmaX1;
    alphaX = -rhoX1*gs/sqrt(gs'*rhoX1*gs);
    bbeta = (g+gX'*(muX1-x)) / sqrt(gs'*rhoX1*gs);
    % deviation from previous bbeta
    deltaBeta = abs(bbeta - index);
    % update beta
    index = bbeta;
    % update x
    x = muX1 + bbeta*sigmaX1.*alphaX;
end

return
end