function [S, q_omp] = OMPDAMAS_DR(Sigma,nbSources,XX, Pmic, k)


% Orthogonal matching pursuit for CMF

%% Initialisation
S=zeros(nbSources,3);

D = [];

% gridded dictionary 
Dom = dictionary(Pmic, XX, k);
normsDR = (sum(abs(Dom).^2, 1).^2 - sum(abs(Dom).^4, 1)).^(1/4);
Dnorm = Dom ./ normsDR;

R = Sigma;

for i = 1:nbSources
%% Computation of the correlations of the residual with each of the atoms of the dictionary,and identification of the maximum

    % gridless beamforming, see below
    [Xbf, d] = bfdr(R,XX, Pmic, k, Dnorm);

%% residual update

    % steering vectors of the identified sources
    D  =[D d];

    % estimation of the powers by CMF NNLS
    q = cmf_nnls_dr(D, Sigma, 0);

    % residual
    R = Sigma - D * diag(q) * D';

    % new position added
    S(i,:)=Xbf;

end
% estimation of the powers
q_omp = cmf_nnls_dr(D, Sigma, 0);

end

% gridless beamforming with diagonal removal
function [Xbf, d] = bfdr(Sigma,XX, Pmic, k, Dnorm)

Sigmadr = Sigma - diag(diag(Sigma));

rho = sum( conj(Dnorm) .*  (Sigmadr * Dnorm), 1);

[~,l_star]=max(rho);

Xbf = fminunc(@(x) obj(x, Pmic, k, Sigma), XX(l_star, :));

d = dictionary(Pmic, Xbf, k);

end

% objective function for beamforming

function f = obj(x, Pmic, k, Sigma)

    d = dictionary(Pmic, x, k);
    D = d*d';
    D = D - diag(diag(D));
    D = D  / norm(D, 'fro');

    f = - abs(trace(D*Sigma));

end
