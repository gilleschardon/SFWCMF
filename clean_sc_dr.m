function [Xest, P, H] = clean_sc_dr(Sigma,Niter,XX, Pmic, k, phi)

% CLEAN-SC
% Sigma SCM
% Niter number of iterations
% XX initialization grid
% Pmic coordinates of the microphones
% k wavenumber
% phi safety factor

Xest = [];
P = [];
H = [];

% normalized steering vectors (for diagonal removal)
Dom = dictionary(Pmic, XX, k);
normsDR = (sum(abs(Dom).^2, 1).^2 - sum(abs(Dom).^4, 1)).^(1/4);
Dnorm = Dom ./ normsDR;

for u = 1:Niter

    % new source by beamforming
    [Xbf, d] = bfdr(Sigma,XX, Pmic, k, Dnorm);

    % source vector
    w = d / sqrt(sum(abs(d).^2)^2 - sum(abs(d).^4));

    % power
    Sigmadr  = Sigma - diag(diag(Sigma));
    p = real(w'*Sigmadr*w);

    Xest = [Xest ; Xbf];
    P = [P p];

    % SC stuff (!!??)
    h = w;
    Sigmadr = Sigma - diag(diag(Sigma));
    for v = 1:100
        HH = diag(diag(h*h'));
        h = (Sigmadr*w/p + HH*w)/sqrt(1 + w'*HH*w);
    end

    H = [H h];
    Sigma = Sigma - phi*p * h * h';
end

end

% see bfdr.m, the steering vectors are here precomputed
function [Xbf, d] = bfdr(Sigma,XX, Pmic, k, Dnorm)

Sigmadr = Sigma - diag(diag(Sigma));
rho = sum( conj(Dnorm) .*  (Sigmadr * Dnorm), 1);
[~,l_star]=max(rho);
Xbf = fminunc(@(x) obj(x, Pmic, k, Sigma), XX(l_star, :));
d = dictionary(Pmic, Xbf, k);

end

function f = obj(x, Pmic, k, Sigma)

    d = dictionary(Pmic, x, k);
    D = d*d';
    D = D - diag(diag(D));
    D = D  / norm(D, 'fro');
    f = - abs(trace(D*Sigma));
end
