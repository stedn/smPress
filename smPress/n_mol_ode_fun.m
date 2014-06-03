

function [ R, birth] = n_mol_ode_fun(q, ti ,R0,Y0,dt)
    %translate input data (time) into reasonable variables
    
    %translate input params into reasonable form
    kon = q(1);
    koff = q(2);
    kph = q(3);
    ks = q(4);
    
    
    [~, val] = ode45(@nmol_ode, ti,[R0 Y0],[], kon, koff, kph, ks);
    R = val(1:length(ti));
    birth = kon*val;
    if(dt~=0)
        [~, val] = ode45(@nmol_ode, 0:(ti(1,2)-ti(1,1)):dt,[R(end) val(end)],[], kon, koff, 0,ks);
        R = [R val];
        birth = [birth kon*val];
    end
end