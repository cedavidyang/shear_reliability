function index = form_GB(MEparam, LD, Rparam)
% return reliabilty index according to partial safty coefficient(PSC)
% PSC: partial safty coefficient
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
L_norm = LD*S/(1.2+1.4*LD);
mu_L = 1.0*L_norm;
sigma_L = 0.25*mu_L;
D_norm = S/(1.2+1.4*LD);
mu_D = 1.05*D_norm;
sigma_D = 0.10*mu_D;

muX = [mu_ME; mu_R; mu_D; mu_L];
sigmaX = [std_ME; std_R; sigma_D; sigma_L];
aEv = sqrt(6)*sigmaX(4)/pi; uEv = -psi(1)*aEv-muX(4);
muX1 = muX; sigmaX1 = sigmaX;
x=muX; normX = eps;
while abs(norm(x)-normX)/normX > 1e-6
    normX = norm(x);
    g = x(1)*x(2)-(x(3)+x(4));
    gX = [x(2); x(1); -1; -1];
    % equivalent normal
    cdfX = [logncdf(x(1), MEparam(1), MEparam(2));
            1-evcdf(-x(4), uEv, aEv)];
    pdfX = [lognpdf(x(1), MEparam(1), MEparam(2));
            evpdf(-x(4), uEv, aEv)];
    nc = norminv(cdfX);
    sigmaX1([1,4]) = normpdf(nc) ./ pdfX;
    muX1([1,4]) = x([1,4]) - nc.*sigmaX1([1,4]);
    gs = gX.*sigmaX1; alphaX = -gs/norm(gs);
    index = (g+gX'*(muX1-x))/norm(gs);
    x = muX1+index*sigmaX1.*alphaX;
end

return
end