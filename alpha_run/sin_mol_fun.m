function [ R ] = sin_mol_fun(q, dat )
    t = dat{1};
%     kon = q(1);
%     koff = q(2);
%     kph = q(3);
%     f = kon + koff + kph;
%     g = sqrt(f^2 - 4*kon*kph);
%     r1 = -f - g;
%     r2 = -f + g;
    r1 = q(1);
    r2 = q(2);
    kph = q(3);
    R0 = q(4);
    A = R0*(kph+r2)/(r2-r1);
    B = -R0*(kph+r1)/(r2-r1);
    R = A*exp(r1*t)+B*exp(r2*t);
end