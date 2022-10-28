function [Xest, Pest, C] = MUSIC_local_1D(Y, nbSources,XX,Pmic, k)

% same as MUSIC_local for a one-dimensional domain

C = Y;

step = XX(2, 2)-XX(1, 2);

[Dom] = dictionary(Pmic, XX, k);

Domnorm = Dom ./ sqrt( sum(abs(Dom).^2, 1));

[V, D] = eigs(C, nbSources, 'largestabs');

ps = sum( abs(V'*Domnorm).^2, 1);

mmax = movmax(ps, 3);

ps(ps ~= mmax) = 0;

idx = find(ps > 0.1);

Xest = [];

tol = 0.001;

if length(idx) > 0
for u = 1:length(idx)

    xopt = fmincon(@(x) obj(x, Pmic, k, V), XX(idx(u), :), [], [], [],[],[],[]);
    zzz = 0;

    for v = 1:size(Xest, 1)
        if norm(xopt - Xest(v, :)) < tol
            zzz = 1;
        end
    end

    if zzz == 0

        Xest = [Xest ; xopt];
    end

end
end



vals = zeros(size(Xest, 1), 1);

for u = 1:size(Xest, 1)
    vals(u) = obj(Xest(u, :), Pmic, k, V);
end

[~, idx] = sort(vals, 'asc');

Xest
vals

nn = min(nbSources, length(idx));

Xest = Xest(idx(1:nn), :);

if numel(Xest) == 0
    Xest = [0,0,0];
    Pest = [0];
    return
end

Dest = dictionary(Pmic, Xest, k);

Dpinv = pinv(Dest);


lambdas = eig(C);
lambdas = sort(abs(lambdas), 'asc');
lambdas(lambdas < 1e-5) = [];
lambda0 = min(lambdas);
lambda0 = mean(lambdas(1:end-nbSources));


PPest = Dpinv * (C - lambda0*eye(size(C))) * Dpinv';

Pest = max(real(diag(PPest)), 0);

end

function f = obj(x, Pmic, k, V)

    d = dictionary(Pmic, x, k);
    d = d / norm(d);

    f = - sum( abs(V'*d).^2, 1);

end
