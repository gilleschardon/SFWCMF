function J = objde(settings, params)

% objective function for differential evolution
Pmic = settings.Pmic;
k = settings.k;
DD = settings.DD;


S1 = dictionary(Pmic, params.X1', k);
S2 = dictionary(Pmic, params.X2', k);
S3 = dictionary(Pmic, params.X3', k);
S4 = dictionary(Pmic, params.X4', k);

C = params.p1 * S1 * S1' + params.p2 * S2 * S2' + params.p3 * S3 * S3' + params.p4 * S4 * S4' + params.noise * eye(size(S1, 1)); 

J = norm(DD - C, 'fro');

end