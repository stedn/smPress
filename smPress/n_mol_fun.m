% Will McFadden (wmcfadden)

function [ R_final ] = n_mol_fun(q, dat )
    %translate input data (time) into reasonable variables
    t = dat{1};
    dt = dat{2};
    R_ret = [];
    birth_ret = [];
    R0 = q(5);
    %Define variables for first time domain solution
    for i=1:length(dt)
        [R, birth] = n_mol_helper(q,t(:,1)',R0,dt(:,i));
        R_ret = [R_ret R(1:length(t(:,1)))];
        birth_ret = [birth_ret birth(1:length(t(:,1))-1)];
        R0 = R(end);
    end
    [R, birth] = n_mol_helper(q,t(:,end)',R0,0);
    R_ret = [R_ret R(1:length(t(:,end)))];
    birth_ret = [birth_ret birth(1:length(t(:,end))-1)];
    R_final = [R_ret birth_ret];
end