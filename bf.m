function [Xbf, d] = bf(Y,XX, Pmic, k)

Dom = dictionary(Pmic, XX, k);
norms = sqrt(sum(abs(Dom.^2), 1));
Domn = Dom ./ norms;

rho = abs(sum(conj(Domn).*(Y*Domn), 1));
[~,l_star]=max(rho);


Xbf = fminunc(@(x) obj(x, Pmic, k, Y), XX(l_star, :));

d = dictionary(Pmic, Xbf, k);

end


function f = obj(x, Pmic, k, Y)

    d = dictionary(Pmic, x, k);
    d = d / norm(d);

    f = - abs(d'*Y*d);
    
end