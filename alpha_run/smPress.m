fitdat = {t};
fitval = R;

r10 =-10*rand;
r20=-10*rand;
kph0 = 10*rand;
rjump = 0.1;
sto = [];
for i = 1:100
    
    param = [r10,r20,kph0];
    lb = [-Inf, -Inf, 0];
    up = [0, 0, Inf];
    [sol,MSE,residual,exitflag,output,lambda,J] = lsqcurvefit(@sin_mol_fun,param,fitdat,fitval,lb,up);

    r1 = sol(1);
    r2 = sol(2);
    k_ph = sol(3);
    R0 = sol(4);

    k_on = r1*r2/k_ph;
    k_off = -(r1+r2)-(k_on+k_ph);
    k_ph;
    sto = [sto; k_on k_off k_ph];
    r10 = r1 + rjump*(rand-0.5);
    r20 = r2 + rjump*(rand-0.5);
    kph0=k_ph+rjump*(rand-0.5);
    rjump = rjump*0.99;
end
covB = inv(J'*J)*MSE;
j_on = [r2/k_ph r1/k_ph -r1*r2/k_ph];
j_off = [1+r2/k_ph 1+r1/k_ph 1 - r1*r2/k_ph];
subcov = covB(1:3,1:3);
var_kon = j_on*subcov*j_on';
var_koff = j_off*subcov*j_off';
var_kph = covB(3,3);
var_R0 = covB(4,4);
plot(nc,xc,'.')
hold on
plot(sin_mol_fun(sol,fitdat),'r')
ylabel('final fit')
finalanswer = full([[k_on; k_off; k_ph; R_0] sqrt([var_kon;var_koff;var_kph;var_R0])])