% Will McFadden (wmcfadden)

function [ R, birth] = n_mol_helper(q, ti ,R0,dt)
    %translate input data (time) into reasonable variables
    
    %translate input params into reasonable form
    r1 = q(1);
    r2 = q(2);
    kph = q(3);
    ks = q(4);
    
    
    %utility variables
    kon = r1*r2/kph;
    koff = -(r1+r2) - (kon+kph);
    %Define variables for first time domain solution
    Ar = (r2+kph)/(r2-r1)*(R0 - ks*r2/kph/(r1+r2+kph));
    Br = (r1+kph)/(r1-r2)*(R0 - ks*r1/kph/(r1+r2+kph));
    
    %Evaluate solution for first time domain 
    R = Ar*exp(r1*ti) + Br*exp(r2*ti) + ks/kph;
    birth = kon*(-Ar*kph/r1*exp(r1*ti(1:end))-Br*kph/r2*exp(r2*ti(1:end))-(r1+r2)*ks/r1/r2 - R(1:end))*(ti(2)-ti(1));
    
    %Evaluate solution at last time point for use in second time domain
    Rt1 = R(end);
    Nt1 = -Ar*kph/r1*exp(r1*ti(end))-Br*kph/r2*exp(r2*ti(end))-(r1+r2)*ks/r1/r2;
    
    %Define variables for second time domain solution
    A = Rt1 + (ks/(r1+r2+kph) + Nt1)*r1*r2/kph/(r1+r2+kph);
    B = -(Nt1+ks/(r1+r2+kph))*r1*r2/kph/(r1+r2+kph);
    
    %Evaluate solution in second time domain
    t2 = 0:ti(1,2)-ti(1,1):dt;
    Rnew = A*exp(-(kon+koff)*t2)+ks*kon*t2/(kon+koff) + B;
    birthnew = kon*(Nt1 + ks*t2(1:end) - Rnew(1:end))*(ti(2)-ti(1));
    R = [R Rnew];
    birth = [birth birthnew];
    
end