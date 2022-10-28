function [Xs, amps, v, nu] = sfw_cmf(Xm, k, Data, Xgrid, lambda, tolpos, tolamp, Niter, LX, UX, Xs, amps, v)


%% SFW for CMF source localization

% Xm microphone positions Mx3
% k wavenumber
% Data covariance matrix, MxM
% X grid initialization grid Nx3
% lambda lambda for l1 regularization
% tolpos tolamp tolerance for source fusion and removal (0 recommended for
% greedy version)
% Niter max number of iterations (in greedy version, = number of sources)
% LX UX bounds of the domain
% Xs amps v, for continuations, solutions for previous lambda, optional

% return
% Xs estimated positions
% amps estimated powers
% v estimated noise variance
% nu stopping criterion

Dgrid = dictionary(Xm, Xgrid, k);
options_nu = optimoptions(@fmincon,'Display', 'off', 'Algorithm','sqp', 'SpecifyObjectiveGradient',true, 'CheckGradient', false, 'OptimalityTolerance', 1e-10);
options_amps = optimoptions(@fmincon,'Display', 'off', 'Algorithm','sqp', 'SpecifyObjectiveGradient',true, 'CheckGradient', false, 'OptimalityTolerance', 1e-10);
options_all = optimoptions(@fmincon,'Display', 'off', 'Algorithm','sqp', 'SpecifyObjectiveGradient',true, 'CheckGradient', false, 'OptimalityTolerance', 1e-10);

if nargin < 11
    Xs = zeros(0, 3);
    amps = zeros(0, 1);
    v = 0;
end

if lambda == 0
    tolnu = 1e-3;
else
    tolnu = 1.01;
end

for u = 1:Niter
    
    fprintf("Iteration %u\t %u sources\t", u, size(Xs, 1))
    
    Dloc = dictionary(Xm, Xs, k);
    C = Dloc*(amps.*Dloc') + v * eye(size(Data));
    
    % residual
    R = Data - C;      
    
    %% Adding a source

    % finding a new source
    [xnew, nu] = maximize_nu(Xm, k, lambda, Xgrid, Dgrid, R, LX, UX, options_nu);

    fprintf("nu = %.8f\n", nu)

    % stop if nu < 1 + eps
    if nu < tolnu
        fprintf("\nnu < tolnu, stopping\n")
        return
    end
    
    % adding the source to the set
    Xs = [Xs ; xnew];

    %% Optimization of the amplitudes
    Dloc = dictionary(Xm, Xs, k);
    [amps, v] = optimize_amplitudes(Dloc, Data, lambda, amps, v, options_amps);  

    %% Joint optimization of the amplitudes and positions 
    [Xs, amps, v] = optimize_all(Xm, k, lambda, Data, Xs, amps, v, LX, UX, options_all);
    
    %% Preparing the next step (cleaning duplicates and 0 power sources)
    [Xs, amps] = clean(Xs, amps, tolamp, tolpos);


end
        %% Max iter. reached
        fprintf("Max iter, stopping\n")
    
    Dloc = dictionary(Xm, Xs, k);
    C = Dloc*(amps.*Dloc');   
end


function [Xnu, nu] = maximize_nu(Xm, k, lambda, Xgrid, Dgrid, R, LX, UX, options)

nugrid = real(sum(conj(Dgrid).*(R*Dgrid), 1));
[~, idx] = max(nugrid);
xnewgrid = Xgrid(idx, :);

nuf = @(X) nux_cov(Xm, k, R, X);
[Xnu, numin] = fmincon(nuf, xnewgrid, [], [], [], [], LX, UX, [], options);

if lambda > 0
nu = -(numin)/lambda;
else
    nu = -numin;
end
end


% removes duplicates and zero power sources
function [Xs, amps] = clean(Xs, amps, tolamp, tolpos)

idxr = amps < tolamp;
amps(idxr) = [];
Xs(idxr, :) = [];

dists = (Xs(:,1) - Xs(:,1)').^2 + (Xs(:,2) - Xs(:,2)').^2 + (Xs(:,3) - Xs(:,3)').^2 + 1000*eye(size(Xs, 1));
[mm, i1] = min(dists);
[m, i2] = min(mm);

i1 = i1(i2);
while m < tolpos    
    xxx = Xs(i1, :);
    Xs([i1 i2], :) = [];
    
    Xs = [Xs ; xxx];
    
    anew = amps(i1, :) + amps(i2, :);
    amps([i1 i2], :) = [];
    
    amps = [amps ; anew];
       
    dists = (Xs(:,1) - Xs(:,1)').^2 + (Xs(:,2) - Xs(:,2)').^2 + (Xs(:,3) - Xs(:,3)').^2 + 1000*eye(size(Xs, 1));
    [mm, i1] = min(dists);
    [m, i2] = min(mm);
    i1 = i1(i2);
end
end

% optimize positions and powers
function [Xs, amps, v] = optimize_all(Xm, k, lambda, Data, Xs, amps, v, LX, UX, options)

    xopt = [Xs(:); amps(:) ; v];

    % bounds
    Ns = length(amps);
    lbounds = [ones(Ns,1)*LX(1) ; ones(Ns, 1)*LX(2); ones(Ns, 1)*LX(3); zeros(Ns+1, 1)];
    ubounds = [ones(Ns,1)*UX(1) ; ones(Ns, 1)*UX(2); ones(Ns, 1)*UX(3); Inf(Ns+1, 1)]; % no upper bounds on amplitudes
    
    
    ZZ = fmincon(@(x) obj_amplitudes_positions(Xm, k, Data, lambda, x), xopt, [], [], [], [], lbounds, ubounds, [], options);

   
    % extaction of the amplitudes and positions
    Xs = reshape(ZZ(1:3*Ns), Ns, 3);
    amps = ZZ(3*Ns+1:end-1);
    v = ZZ(end);
    
end

% optimize powers
function [amps, v] = optimize_amplitudes(Dloc, Data, lambda, amps, v, options)

    [ampsv] = fmincon(@(x) obj_amplitudes(Dloc, Data, lambda, x), [amps ; eps ; v], [], [], [], [], zeros(size(amps, 1)+2,1),[], [], options);
    amps = ampsv(1:end-1);
    v = ampsv(end);
     
end


%% objective functions (include gradients)

function [J, Jgrad] = obj_amplitudes_positions(Xm, k, Data, lambda, xx)

Ns = (length(xx)-1)/4;

Xs = reshape(xx(1:3*Ns), Ns, 3);
x = xx(3*Ns+1:end-1);
v = xx(end);

if nargout > 1

[D, Dx, Dy, Dz] =  dictionary(Xm, Xs, k);

C = D*(x.*D') + v * eye(size(Data));

Delta = C - Data;

J = 1/2 * norm(Delta, 'fro').^2 + lambda * sum(x);

JD = Delta*D;

    Jgrada = ones(size(x)) * lambda;
    Jgrada = Jgrada + real(sum(JD .* conj(D), 1))'; 

    Jgradx = 2 * x .* real(sum(JD .* conj(Dx), 1))';
    Jgrady = 2 * x .* real(sum(JD .* conj(Dy), 1))';
    Jgradz = 2 * x .* real(sum(JD .* conj(Dz), 1))';

    Jgradv = real(trace(Delta));

    Jgrad = [Jgradx ; Jgrady ; Jgradz ; Jgrada ; Jgradv];

else
    
D =  dictionary(Xm, Xs, k);
C = D*(x.*D') + v * eye(size(Data));
Z = C - Data;

J = 1/2 * norm(Z, 'fro').^2 + lambda * sum(x); 
end
end


function [J, Jgrad] = obj_amplitudes(D, Data, lambda, x)

len = length(x);

C = D*(x(1:end-1).*D') + eye(size(Data)) * x(end);

Delta = C - Data;

J = 1/2 * norm(Delta, 'fro').^2 + lambda * sum(x(1:end-1));

Jgrad = zeros(len, 1);

if nargout > 1
    Jgrad(1:end-1) = ones(size(x(1:end-1))) * lambda;
    Jgrad(1:end-1) = Jgrad(1:end-1) + real(sum( (Delta*D) .* conj(D), 1))'; 
    Jgrad(end) = real(trace(Delta));
end

end

function [nu, nugrad] = nux_cov(Xm, k, Delta, Xnu)


if nargout > 1
    [d, gx, gy, gz] = dictionary(Xm, Xnu, k);
    d = d(:);
    
    dD = d'*(Delta);
    nu = -real(dD*d);
    
    ggx = - 2 * real(dD*gx);
    ggy = - 2 * real(dD*gy);
    ggz = - 2 * real(dD*gz);
  
    nugrad = [ggx ggy ggz];
    

else
    d = dictionary(Xm, Xnu, k);

    nu = -real(d'*(Delta)*d);

end
end
