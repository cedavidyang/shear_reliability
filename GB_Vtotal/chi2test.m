function [h p] = chi2test(x, nbins, alpha, parameter, nparams, distribution)
% chi-square test of sample x
% x: sample colomn vector
% nbins: number of bins
% alpha: level of significance
% parameter: parameters of assumed distribution
% nparams: number of estimated parameters, it may or may not be the same as
% the dimension of parameter
% distribution: name of the distribution

n = length(x);
% adjust nbins
while 1
    cdf = (1:nbins)/nbins;
    E = zeros(nbins, 1);
    for i=1:nbins
        if i~=nbins
            E(i) = 1/nbins*n;
        else
            E(i) = n-(nbins-1)*E(i-1);
        end
    end
    O = zeros(nbins,1);
    switch distribution
        case 'Gamma'
            O(1) = length( x( x>gaminv(0, parameter(1), parameter(2)) & x<gaminv(cdf(1), parameter(1), parameter(2)) ) );
            for i=2:nbins
                O(i) = length( x( x>gaminv(cdf(i-1), parameter(1), parameter(2)) & x<gaminv(cdf(i), parameter(1), parameter(2)) ) );
            end
        case 'Lognormal'
            % note: parameters are MEAN and DEVIATION of associated normal
            % distribution
            O(1) = length( x( x>logninv(0, parameter(1), parameter(2)) & x<logninv(cdf(1), parameter(1), parameter(2)) ) );
            for i=2:nbins
                O(i) = length( x( x>logninv(cdf(i-1), parameter(1), parameter(2)) & x<logninv(cdf(i), parameter(1), parameter(2)) ) );
            end
        case 'Weibull'
            O(1) = length( x( x>wblinv(0, parameter(1), parameter(2)) & x<wblinv(cdf(1), parameter(1), parameter(2)) ) );
            for i=2:nbins
                O(i) = length( x( x>wblinv(cdf(i-1), parameter(1), parameter(2)) & x<wblinv(cdf(i), parameter(1), parameter(2)) ) );
            end
        case 'Normal'
            O(1) = length( x( x>norminv(0, parameter(1), parameter(2)) & x<norminv(cdf(1), parameter(1), parameter(2)) ) );
            for i=2:nbins
                O(i) = length( x( x>norminv(cdf(i-1), parameter(1), parameter(2)) & x<norminv(cdf(i), parameter(1), parameter(2)) ) );
            end
        otherwise
            disp('unknown distribution')
    end
    if ~isempty( O(O<5) )
        nbins = nbins-1;
        if nbins < 4+nparams
            break;
        end
    else
        break;
    end
end

if nbins == 3+nparams
    disp('sample is too small')
else
    bins = 0:(nbins-1);
    [h,p,st] = chi2gof(bins,'ctrs',bins,'frequency',O,'expected',E, 'nparams',nparams, 'alpha', alpha);
end

return
end
