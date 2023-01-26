function [Xest, P, H, Pcsc] = hr_clean_sc(Y,Niter,XX, Pmic, k, phi, mu)

Xest = [];
P = [];
H = [];
Pcsc = [];
G = [];
Y0 = Y;
W = [];
Dom = dictionary(Pmic, XX, k);

normsDR = (sum(abs(Dom).^2, 1).^2 - sum(abs(Dom).^4, 1)).^(1/4);

Dnorm = Dom ./ normsDR;



for u = 1:Niter
    [Xbf, d] = bf(Y,XX, Pmic, k);

    w = d / norm(d)^2;
       
    

    p = real(w'*Y*w);
    
    
    
    Xest = [Xest ; Xbf];
    P = [P p];

    h = w;

    h = Y * w / (w'*Y*w);
    
    H = [H h];
    G = [G d];
    Y = Y - phi*p * (h * h');


end

U = G ./ sum(abs(G).^2, 1);
NU = zeros(size(U));
NiterHR = 10;

for v = 1:NiterHR
    for u = 1:Niter
        if mu == 0
            nu = fminunc(@(XRI) funHR(G, u, XRI), [real(U(:, u)) ; imag(U(:, u))]);
        else
            nu = fmincon(@(XRI) funHR(G, u, XRI), [real(U(:, u)) ; imag(U(:, u))], [], [], [], [], [], [], @(XRI) nonlcon(XRI, G(:, u), mu));
        end
        
        l = length(nu)/2;
        NU(:, u) = nu(1:l) + 1i * nu(l+1:end);
    end
    U = NU;
    for u = 1:Niter
        H(:, u) = Y0 * U(:, u) / (U(:, u)' * Y0 * U(:, u));
        Cu = H(:, u) * H(:, u)';
        [Xbf, d] = bf(Cu,XX, Pmic, k);
        G(:, u) = d;
        Xest(u, :) = Xbf;
        P(u) = real(U(:, u)'*Y0*U(:, u)) * abs(d'*H(:, u)/norm(d)^2)^2;
        
        
    end
    
    
end
            


end

function [c, ceq] = nonlcon(u, g, mu)
    ceq = [];
    l = length(u)/2;
    u = u(1:l) + 1i*u(l+1:end);
    c = (mu - abs(g'*u)^2);

end

function obj = funHR(G, idx, u)
    l = length(u)/2;
    u = u(1:l) + 1i*u(l+1:end);
    vecs = G * diag(G'*u);
    vecs(:, idx) = 0;
    N = sum(abs(sum(vecs, 2)).^2);
    
    obj = N / ( abs(G(:, idx)'*u)^2 * norm(G(:, idx))^2);
    
end