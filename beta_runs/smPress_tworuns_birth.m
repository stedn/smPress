

% produce the true values up front
kont = 0.005;
kofft = 0.1;
kpht = 0.1;
kst = 1.0;
R0t = 200;
r1t = -(kont+kofft+kpht)/2 - sqrt(((kont+kofft+kpht))^2-4*kont*kpht)/2;
r2t = -(kont+kofft+kpht)/2 + sqrt(((kont+kofft+kpht))^2-4*kont*kpht)/2;

%produce a couple of time domains starting from 0
t1 = 0:100;
t2 = 0:150;

%produce a wait time
dt = 180;

%initialize the fake data
%formatting
pinit = [r1t,r2t,kpht,kst,R0t];
datinit = {t1,t2,dt};
%call function
Rinit = dub_mol_cyto_fun(pinit,datinit);
%split up
R1 = Rinit(1:length(t1));
R2 = Rinit(length(t1)+1:length(t1)+length(t2));
birth1 = Rinit(length(t1)+length(t2)+1:2*length(t1)-1+length(t2));
birth2 = Rinit(2*length(t1)+length(t2):end);

noise = 5;
R1 = R1+noise*(rand(size(R1))-0.5);
R2 = R2+noise*(rand(size(R2))-0.5);
birth1 = birth1+0.1*noise*(rand(size(birth1))-0.5);
birth2 = birth2+0.1*noise*(rand(size(birth2))-0.5);

%formatting for call to fitter
fitdat = {t1, t2, dt};
fitval = [R1 R2 birth1 birth2];

%randomish initial guesses for fit params
r10 =0.5*r1t+r1t*rand;
r20=0.5*r2t+r2t*rand;
kph0 = 0.5*kpht+kpht*rand;
ks0 = 0.5*kst+kst*rand;
R00=0.5*R0t+R0t*rand;

%for repeated stability of fit test
rjump = 1;
sto = [];
stosol = [];
for i = 1:10
    
    %call fitter
    param = [r10,r20,kph0,ks0,R00];
    lb = [-Inf, -Inf, 0, 0, R00*0.5];
    up = [0, 0, Inf, Inf, R00*2];
    [sol,MSE,residual,exitflag,output,lambda,J] = lsqcurvefit(@dub_mol_birth_fun,param,fitdat,fitval,lb,up);
    
    %rename fit params
    r1 = sol(1);
    r2 = sol(2);
    k_ph = sol(3);
    k_s = sol(4);
    R_0 = sol(5);

    %evaluate reality params from r1 r2
    k_on = r1*r2/k_ph;
    k_off = -(r1+r2)-(k_on+k_ph);
    %store info
    sto = [sto; k_on k_off k_ph k_s R_0];
    
    %random initial guesses from output
    r10 = r1 + r1*rjump*(rand-0.5);
    r20 = r2 + r2*rjump*(rand-0.5);
    kph0=k_ph+k_ph*rjump*(rand-0.5);
    ks0=k_s+k_s*rjump*(rand-0.5);
    R00=R_0+R_0*rjump*(rand-0.5);
    
    %decrease annealing temp
    rjump = rjump*0.99;
end
covB = inv(J'*J)*MSE;
j_on = [r2/k_ph r1/k_ph -r1*r2/k_ph];
j_off = [1+r2/k_ph 1+r1/k_ph 1 - r1*r2/k_ph];
subcov = covB(1:3,1:3);
var_kon = j_on*subcov*j_on';
var_koff = j_off*subcov*j_off';
var_kph = covB(3,3);
var_ks = covB(4,4);
var_R0 = covB(5,5);
    
% [beta,R,J,covB,MSE] = nlinfit(fitdatbeta,fitval,@dub_mol_fun_beta,sol);
    
%plot everything
% subplot(2,3,1)
% plot(sto(:,1))
% ylabel('k_{on}')
% subplot(2,3,2)
% plot(sto(:,2))
% ylabel('k_{off}')
% subplot(2,3,3)
% plot(sto(:,3))
% ylabel('k_{ph}')
% subplot(2,3,4)
% plot(sto(:,4))
% ylabel('k_{s}')
% subplot(2,3,5)
% plot(sto(:,5))
% ylabel('R_{0}')
fsol = dub_mol_cyto_fun(sol,fitdat);
% subplot(2,1,1)
plot(t1,R1,'.',t1(end)+dt+t2,R2,'.')
hold on
plot([t1 (t2+t1(end))],fsol(1:length(t1)+length(t2)),'r')
ylabel('final fit to R')
% subplot(2,1,2)
% plot(t1(2:end),birth1,'.',t1(end)+t2(2:end),birth2,'.')
% hold on
% plot([t1(2:end) (t2(2:end)+t1(end))],fsol(length(t1)+length(t2)+1:end),'r')
% ylabel('final fit to R')
finalanswer = full([[k_on; k_off; k_ph; k_s; R_0] sqrt([var_kon;var_koff;var_kph;var_ks;var_R0])])