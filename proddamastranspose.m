function x = proddamastranspose(D, Data)

% fast product by Dtilde transpose
% see
% G. Chardon, J. Picheral, F. Ollivier
% Theoretical analysis of the DAMAS algorithm and efficient implementation of the Covariance Matrix Fitting method for large-scale problems
% Journal of Sound and Vibration, 2021
% doi:10.1016/j.jsv.2021.116208
x = real(sum((D'*Data).*D.', 2));

end