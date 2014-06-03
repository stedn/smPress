function [ R ] = dub_mol_birth_fun(q, dat )
    %translate input data (time) into reasonable variables
    t1 = dat{1};
    t2 = dat{2};
    dt = dat{3};

    %translate input params into reasonable form
    r1 = q(1);
    r2 = q(2);
    kph = q(3);
    ks = q(4);
    R0 = q(5);
    
    %utility variables
    kon = r1*r2/kph;
    koff = -(r1+r2) - (kon+kph);
    
    %Define variables for first time domain solution
    Ar = (r2+kph)/(r2-r1)*(R0 - ks*r2/kph/(r1+r2+kph));
    Br = (r1+kph)/(r1-r2)*(R0 - ks*r1/kph/(r1+r2+kph));
    
    %Evaluate solution for first time domain 
    R1_ret = Ar*exp(r1*t1) + Br*exp(r2*t1) + ks/kph;
    birth1_ret = kon*(-Ar*kph/r1*exp(r1*t1(2:end))-Br*kph/r2*exp(r2*t1(2:end))-(r1+r2)*ks/r1/r2 - R1_ret(2:end));
%     birth1_ret = kon*(-Ar*(1+kph/r1)*exp(r1*t1)-Br*(1+kph/r2)*exp(r2*t1)-ks/kph-(r1+r2)/(r1*r2)*ks);
    %Evaluate solution at last time point for use in second time domain
    Rt1 = R1_ret(end);
    Nt1 = -Ar*kph/r1*exp(r1*t1(end))-Br*kph/r2*exp(r2*t1(end))-(r1+r2)*ks/r1/r2;
    
    %Define variables for second time domain solution
    A = Rt1 + (ks/(r1+r2+kph) + Nt1)*r1*r2/kph/(r1+r2+kph);
    B = -(Nt1+ks/(r1+r2+kph))*r1*r2/kph/(r1+r2+kph);
    
    %Evaluate solution at last time point for use in third time domain
    Rt2 = A*exp(-(kon+koff)*dt)+ks*kon*dt/(kon+koff) + B;
    
    %Define variables for third time domain solution
    Cr = (r2+kph)/(r2-r1)*(Rt2 - ks*r2/kph/(r1+r2+kph));
    Dr = (r1+kph)/(r1-r2)*(Rt2 - ks*r1/kph/(r1+r2+kph));
    
    %Evaluate solution for third time domain 
    R2_ret = Cr*exp(r1*t2) + Dr*exp(r2*t2) + ks/kph;
    birth2_ret = kon*(-Cr*kph/r1*exp(r1*t2(2:end)) - Dr*kph/r2*exp(r2*t2(2:end)) - (r1+r2)*ks/r1/r2 - R2_ret(2:end));
%     birth2_ret = kon*(-Cr*(1+kph/r1)*exp(r1*t2)-Dr*(1+kph/r2)*exp(r2*t2)-ks/kph-(r1+r2)/(r1*r2)*ks);
    
    %return concatenated solution for evaluation by least squares fit
    R = [R1_ret R2_ret birth1_ret birth2_ret];
end