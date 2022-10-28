function y = proddamas_support(D, x, S)

% fast product by Dtilde transpose D (product by the large DAMAS matrix)
% see
% G. Chardon, J. Picheral, F. Ollivier
% Theoretical analysis of the DAMAS algorithm and efficient implementation of the Covariance Matrix Fitting method for large-scale problems
% Journal of Sound and Vibration, 2021
% doi:10.1016/j.jsv.2021.116208

z = D(:, S)*(x(S).* D(:, S)');
y = real(sum(conj(D) .* (z*D), 1).');


end