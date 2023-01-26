function [Xest, P, H, Pcsc] = hr_clean_sc_dr(Y,Niter,XX, Pmic, k, phi, mu)

Xest = [];
P = [];
H = [];
Pcsc = [];
G = [];
Y0 = Y;
Y0dr = Y0 - diag(diag(Y0));
Dom = dictionary(Pmic, XX, k);

normsDR = (sum(abs(Dom).^2, 1).^2 - sum(abs(Dom).^4, 1)).^(1/4);

Dnorm = Dom ./ normsDR;



for u = 1:Niter
    [Xbf, d] = bf(Y,XX, Pmic, k);

    w = d / sqrt(sum(abs(d).^2)^2 - sum(abs(d).^4));

    
    Ydr  = Y - diag(diag(Y));
    p = real(w'*Ydr*w);
    
    
    

    Xest = [Xest ; Xbf];
    P = [P p];

    
    h = d;
    Ydr = Y - diag(diag(Y));
    for v = 1:100
        HH = diag(diag(h*h'));
        h = (Ydr*w/p + HH*w)/sqrt(1 + w'*HH*w);
    end
    
    H = [H h];
    G = [G d];
    Y = Y - phi*p * (h * h');


end

U = G ./ sqrt(sum(abs(G).^2, 1).^2 - sum(abs(G).^4, 1));
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
        H(:, u) = G(:, u);
        
        
        for v = 1:100
            HH = diag(diag(H(:, u)*H(:, u)'));
            h = (Y0dr*U(:, u)/(U(:, u)'*Y0dr*U(:, u)) + HH*U(:, u))/sqrt(1 + U(:, u)'*HH*U(:, u));
        end
        H(:, u) = h;


        Cu = H(:, u) * H(:, u)';
        [Xbf, d] = bfdr(Cu,XX, Pmic, k, Dnorm);
        
        Xest(u, :) = Xbf;
        G(:, u) = d;

        
        
        
        
        HH = H(:, u)*H(:, u)';
        HH = HH - diag(diag(HH));
        w = d / sqrt(sum(abs(d).^2)^2 - sum(abs(d).^4));

        P(u) = real(U(:, u)'*Y0dr*U(:, u)) * abs(w'*HH*w);
        
        
    end
    
    
end
            


end


function [c, ceq] = nonlcon(u, g, mu)
    ceq = [];
    l = length(u)/2;
    u = u(1:l) + 1i*u(l+1:end);
    c = (mu - abs(g'*u)^2);

end

function [Xbf, d] = bfdr(Y,XX, Pmic, k, Dnorm)

Ydr = Y - diag(diag(Y));
rho = real(sum( conj(Dnorm) .*  (Ydr * Dnorm), 1));
[~,l_star]=max(rho);
Xbf = fminunc(@(x) obj(x, Pmic, k, Y), XX(l_star, :));


d = dictionary(Pmic, Xbf, k);

end

function f = obj(x, Pmic, k, Y)

    d = dictionary(Pmic, x, k);
    D = d*d';
    D = D - diag(diag(D));
    D = D  / norm(D, 'fro');
    f = - abs(trace(D*Y));

end



function obj = funHR(G, idx, u)
    l = length(u)/2;
    u = u(1:l) + 1i*u(l+1:end);
    normsG = sum(abs(G).^2, 1);
    scgu = abs(u' * G).^2;
    vecs = scgu .* normsG;
    vu = vecs(idx);
    vecs(idx) = 0;
    N = sum(abs(sum(vecs, 2)).^2);
    
    obj = N / vu;
    
end