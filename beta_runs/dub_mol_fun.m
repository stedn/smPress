function [ R ] = dub_mol_fun(q, dat )
    %translate input data (time) into reasonable variables
    t1 = dat{1};
    t2 = dat{2};
    dt = dat{3};

    %translate input params into reasonable form
    r1 = q(1);
    r2 = q(2);
    kph = q(3);
    R0 = q(4);
    
    %Define variables for first time domain solution
    Ar = (r2+kph)/(r2-r1)*(R0 );
    Br = (r1+kph)/(r1-r2)*(R0 );
    
    %Evaluate solution for first time domain 
    R1_ret = Ar*exp(r1*t1) + Br*exp(r2*t1) ;
    
    %Evaluate solution at last time point for use in second time domain
    Rt1 = R1_ret(end);
    Nt1 = -Ar*kph/r1*exp(r1*t1(end))-Br*kph/r2*exp(r2*t1(end));
    
    %Define variables for second time domain solution
    A = Rt1 + ( Nt1)*r1*r2/kph/(r1+r2+kph);
    B = -(Nt1)*r1*r2/kph/(r1+r2+kph);
    kon = r1*r2/kph;
    koff = -(r1+r2) - (kon+kph);
    
    %Evaluate solution at last time point for use in third time domain
    Rt2 = A*exp(-(kon+koff)*dt) + B;
    
    %Define variables for third time domain solution
    Cr = (r2+kph)/(r2-r1)*(Rt2 );
    Dr = (r1+kph)/(r1-r2)*(Rt2 );
    
    %Evaluate solution for third time domain 
    R2_ret = Cr*exp(r1*t2) + Dr*exp(r2*t2) ;
    
    %return concatenated solution for evaluation by least squares fit
    R = [R1_ret R2_ret];
end