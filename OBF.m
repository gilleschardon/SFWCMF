function [XS, P] = OBF(Sigma,nbSources,XX, Pmic, k)

% orthogonal beamforming

% eigendecomposition
[U, S, V] = svd(Sigma);

P = zeros(nbSources, 1);
XS = zeros(nbSources, 3);

for i = 1:nbSources

    % gridless beamforming
    [Xloc, d] = bf(U(:, i)*U(:, i)',XX, Pmic, k);

    % power
    rho = abs((d'/norm(d)*U(:, i) * sqrt(S(i,i))).^2);

    XS(i, :) = Xloc;
    P(i) = rho/norm(d)^2;

end


end
