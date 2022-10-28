function [Xbf, d] = bfdr(Sigma,XX, Pmic, k)

% beamforming with diagonal removal, with local optimization
% Sigma sample covariance matrix
% XX coordinates of the init grid
% Pmic coordinate of the microphones
% k wavenumber

% Xbf estimated position
% d steering vector at Xbf

% normalized steering vectors

Dom = dictionary(Pmic, XX, k);
normsDR = (sum(abs(Dom).^2, 1).^2 - sum(abs(Dom).^4, 1)).^(1/4);
Dnorm = Dom ./ normsDR;

% diagonal removal
Sigmadr = Sigma - diag(diag(Sigma));

% beamforming with diagonal removal (note that the diagonal is removed in the data only, however the steering vectors have the appropriate normalization)
rho = sum( conj(Dnorm) .*  (Sigmadr * Dnorm), 1);

% grid point for the initialization
[~,l_star]=max(rho);

% local optimization
Xbf = fminunc(@(x) obj(x, Pmic, k, Sigma), XX(l_star, :));

% return the source vector
d = dictionary(Pmic, Xbf, k);

end

% beamforming objective function
function f = obj(x, Pmic, k, Sigma)

    d = dictionary(Pmic, x, k);
    D = d*d';
    D = D - diag(diag(D));
    D = D  / norm(D, 'fro');

    f = - abs(trace(D*Sigma));

end
