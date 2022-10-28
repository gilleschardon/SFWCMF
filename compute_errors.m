function [epos, eamps] = compute_errors(Xest, Xtrue, Aest, Atrue)

% compute position and amplitude errors by matching sources
% does not check if all sources have been found

epos = 0;
eamps = 0;

if numel(Xest) > 0

dists = (Xest(:,1) - Xtrue(:,1)').^2 + (Xest(:,2) - Xtrue(:,2)').^2 + (Xest(:,3) - Xtrue(:,3)').^2;

M = matchpairs(dists, 1e10);



for u = 1:size(M, 1)
    epos = epos + norm(Xest(M(u, 1), :) - Xtrue(M(u, 2), :))^2;
    eamps = eamps + norm(abs(Aest(M(u, 1))) - abs(Atrue(M(u, 2))), 'fro')^2;
end

end

end
