function index = sorm_HK(PSC, MEparam, LD, Rparam)
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
norm_R = Rparam(3);

S = norm_R*PSC;
L_norm = LD*S/(1.4+1.6*LD);
% L_norm = LD*S/(1.4+1.7*LD);
% L_norm = LD*S/(1.2+1.4*LD);
mu_L = 1.2*L_norm;
std_L = 0.25*mu_L;
D_norm = S/(1.4+1.6*LD);
% D_norm = S/(1.4+1.7*LD);
% D_norm = S/(1.2+1.4*LD);
mu_D = 1.05*D_norm;
std_D = 0.10*mu_D;

muX = [mu_ME; mu_R; mu_D; mu_L];
sigmaX = [std_ME; std_R; std_D; std_L];
aEv = sqrt(6)*sigmaX(4)/pi; uEv = -psi(1)*aEv-muX(4);
muX1 = muX; sigmaX1 = sigmaX;
x=muX; normX = eps;
cdfX = [logncdf(x(1), MEparam(1), MEparam(2)); normcdf( x(2), mu_R, std_R );...
        normcdf(x(3), mu_D, std_D); 1-evcdf(-x(4), uEv, aEv)];
y = norminv(cdfX);    
while abs(norm(x)-normX)/normX > 1e-6
    normX = norm(x);
    g = x(1)*x(2)-(x(3)+x(4));
    gX = [x(2); x(1); -1; -1];
    pdfX = [lognpdf(x(1), MEparam(1), MEparam(2)); normpdf( x(2), mu_R, std_R );...
        normpdf(x(3), mu_D, std_D); 1-evpdf(-x(4), uEv, aEv)];
    xY = normpdf(y) ./ pdfX; gY = gX.*xY;
    alphaY = -gY/norm(gY);
    bbeta = (g-gY'*y) / norm(gY);
    y = bbeta*alphaY;
    cdfY = normcdf(y);
    x = [logninv(cdfY(1), MEparam(1), MEparam(2)); norminv( cdfY(2), mu_R, std_R );...
        norminv(cdfY(3), mu_D, std_D); -evinv(1-cdfY(4), uEv, aEv)];    
end %FORM to SORM
n=length(muX);
h=[null(alphaY'), alphaY];
gXX = [0, 1, 0, 0; 1, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0];
mLn = MEparam(1); sLn = MEparam(2);
fX = [(mLn-sLn^2-log(x(1)))/sLn^2/x(1); (mu_R-x(2))/(std_R^2);...
      (mu_D-x(3))/(std_D^2); (exp((-uEv-x(4))/aEv)-1)/aEv].*pdfX;
gYY = xY*xY' .*gXX - diag( gX.*(y.*xY+xY.^2.*fX./pdfX) );
q = -gYY / norm(gY);
qt = h'*q*h; qt(:, n) = []; qt(n, :) = [];
pFQ = normcdf(-bbeta) / sqrt(det(eye(n-1)-bbeta*qt));
index = -norminv( pFQ, 0, 1);
% ME_ptile = 1;
% muX = [mu_R; mu_D; mu_L];
% sigmaX = [std_R; std_D; std_L];
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

return
end