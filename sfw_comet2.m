function [Xs, amps, v, nu, Pest] = sfw_comet2(Xm, k, Data, Xgrid, tolpos, tolamp, Niter, LX, UX)

%% SFW for COMET2

% Xm microphone positions Mx3
% k wavenumber
% Data covariance matrix, MxM
% Xgrid initialization grid Nx3
% tolpos tolamp tolerance for source fusion and removal (0 recommended for
% greedy version)
% Niter max number of iterations (in greedy version, = number of sources)
% LX UX bounds of the domain

% return
% Xs estimated positions
% amps estimated powers
% v estimated noise variance
% nu stopping criterion
% Pest power reestimated by least-squares

Dgrid = dictionary(Xm, Xgrid, k);
options_nu = optimoptions(@fmincon,'Display', 'off', 'Algorithm','sqp', 'OptimalityTolerance', 1e-12);
options_amps = optimoptions(@fmincon,'Display', 'off', 'Algorithm','sqp', 'SpecifyObjectiveGradient',false, 'OptimalityTolerance', 1e-12);
options_all = optimoptions(@fmincon,'Display', 'off', 'Algorithm','sqp',  'OptimalityTolerance', 1e-12);

Xs = zeros(0, 3);
amps = zeros(0, 1);
v = 1;

lbv = 0.01;

for u = 1:Niter
    
    fprintf("Iteration %u\t %u sources\t", u, size(Xs, 1))
    
    Dloc = dictionary(Xm, Xs, k);
    C = Dloc*(amps.*Dloc') + (lbv+v) * eye(size(Data));
    
    Cinv = eye(size(Data))/(lbv + v) - Dloc * inv(diag(1./amps) + Dloc'*Dloc/(lbv + v)) * Dloc' / (lbv+v) ^2;
    
    %% Adding a source

    [xnew, nu] = maximize_nu(Xm, k, Xgrid, Dgrid, Data, C, Cinv, LX, UX, options_nu);
    
    fprintf("nu = %.4f\n", nu)


    
    % adding the source to the set
    Xs = [Xs ; xnew];

    %% Optimization of the amplitudes
    Dloc = dictionary(Xm, Xs, k);

    
    [amps, v] = optimize_amplitudes(Dloc, Data, amps, v, options_amps, lbv);  

    %% Joint optimization of the amplitudes and positions 
    [Xs, amps, v] = optimize_all(Xm, k, Data, Xs, amps, v, LX, UX, options_all, lbv);
    
    %% Preparing the next step
    [Xs, amps] = clean(Xs, amps, tolamp, tolpos);


end

D = dictionary(Xm, Xs, k);
Dpinv = pinv(D);



lambdas = eig(Data);
lambdas = sort(abs(lambdas), 'asc');
lambda0 = mean(lambdas(1:end-Niter));

PPest = Dpinv * (Data - lambda0*eye(size(Data))) * Dpinv';



Pest = max(real(diag(PPest)), 0);
        %% Max iter. reached
        fprintf("Max iter, stopping\n")
    
end


function [Xnu, nu] = maximize_nu(Xm, k, Xgrid, Dgrid, Data, C, Cinv, LX, UX, options)

RR = Cinv*Data*Data*Cinv;

nugrid = real(sum(conj(Dgrid).*(RR*Dgrid), 1)) - (sum(abs(Dgrid).^2, 1));

[~, idx] = max(nugrid);
xnewgrid = Xgrid(idx, :);



nuf = @(X) nux_cov(Xm, k, RR, X);
[Xnu, numin] = fmincon(nuf, xnewgrid, [], [], [], [], LX, UX, [], options);
nu = -numin;
%Xnu
end


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

function [Xs, amps, v] = optimize_all(Xm, k,Data, Xs, amps, v, LX, UX, options, lbv)

    xopt = [Xs(:); amps(:) ; v];

    % bounds
    Ns = length(amps);
    lbounds = [ones(Ns,1)*LX(1) ; ones(Ns, 1)*LX(2); ones(Ns, 1)*LX(3); zeros(Ns+1, 1)];
    ubounds = [ones(Ns,1)*UX(1) ; ones(Ns, 1)*UX(2); ones(Ns, 1)*UX(3); Inf(Ns+1, 1)]; % no upper bounds on amplitudes
    
    
    ZZ = fmincon(@(x) obj_amplitudes_positions(Xm, k, Data, x, lbv), xopt, [], [], [], [], lbounds, ubounds, [], options);

   
    % extaction of the amplitudes and positions
    Xs = reshape(ZZ(1:3*Ns), Ns, 3);
    amps = ZZ(3*Ns+1:end-1);
    v = ZZ(end);
    
end

function [amps, v] = optimize_amplitudes(Dloc, Data, amps, v, options, lbv)

    [ampsv] = fmincon(@(x) obj_amplitudes(Dloc, Data, x, lbv), [amps ; eps ; v], [], [], [], [], zeros(size(amps, 1)+2,1),[], [], options);
    amps = ampsv(1:end-1);
    v = ampsv(end);
     
end

function [J] = obj_amplitudes_positions(Xm, k, Data, xx, lbv)

Ns = (length(xx)-1)/4;

Xs = reshape(xx(1:3*Ns), Ns, 3);
x = xx(3*Ns+1:end-1);
v = xx(end);

D =  dictionary(Xm, Xs, k);
C = D*(x.*D') + (lbv+v) * eye(size(Data));

sup = (x >= 10*eps);

Z = diag(1./x(sup)) + D(:, sup)'*D(:, sup)/(lbv + v);


Cinv = eye(size(Data))/(lbv + v) - D(:, sup) * inv(Z) * D(:, sup)' / (lbv+v) ^2;

J = real(trace(Data*Cinv*Data) + trace(C));


end



function [J, Jgrad] = obj_amplitudes(D, Data, x, lbv)

len = length(x);

C = D*(x(1:end-1).*D') + eye(size(Data)) * (lbv+x(end));

a = x(1:end-1);
v = x(end);

sup = (a > 10*eps);


Z = diag(1./a(sup)) + D(:, sup)'*D(:, sup)/(lbv + x(end));


Cinv = eye(size(Data))/(lbv + x(end)) - D(:, sup) * inv(Z) * D(:, sup)' / (lbv+x(end)) ^2;


J = real(trace(Data*Cinv*Data) + trace(C));


end

function [nu] = nux_cov(Xm, k, RR, Xnu)

    d = dictionary(Xm, Xnu, k);
    
    

    nu = - real(d'*(RR)*d) + norm(d)^2;

end

