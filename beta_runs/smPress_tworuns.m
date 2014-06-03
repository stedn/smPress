

% produce the true values up front
kont = 0.01;
kofft = 0.08;
kpht = 0.1;
kst = 0.0;
R0t = 400;
r1t = -(kont+kofft+kpht)/2 - sqrt(((kont+kofft+kpht))^2-4*kont*kpht)/2;
r2t = -(kont+kofft+kpht)/2 + sqrt(((kont+kofft+kpht))^2-4*kont*kpht)/2;

%produce a couple of time domains starting from 0
% t1 = 0:0.5:50;
% t2 = 0:0.5:50;

%produce a wait time
dt = 30;

%initialize the fake data
%formatting
pinit = [r1t,r2t,kpht,kst,R0t];
datinit = {t1,t2,dt};
%call function
Rinit = dub_mol_fun(pinit,datinit);
%split up
% R1 = Rinit(1:length(t1));
% R2 = Rinit(length(t1)+1:end);

% noiser = 0.25;
% R1 = R1.*(1+noiser*(rand(size(R1))-0.5));
% R2 = R2.*(1+noiser*(rand(size(R2))-0.5));

 
%formatting for call to fitter
fitdat = {t1 t2 dt};
fitval = [R1 R2];

%randomish initial guesses for fit params
r10 =r1t;
r20=r2t;
kph0 = kpht;
ks0 = 0.1;
R00=R1(1);

%for repeated stability of fit test
rjump = 0.1;
sto = [];
for i = 1:100
    
    %call fitter
    param = [r10,r20,kph0,R00];
    lb = [-Inf, -Inf, 0,  R00*0.5];
    up = [0, 0, Inf,  R00*2];
    [sol res] = lsqcurvefit(@dub_mol_fun,param,fitdat,fitval,lb,up);

    %rename fit params
    r1 = sol(1);
    r2 = sol(2);
    k_ph = sol(3);
    R_0 = sol(4);

    %evaluate reality params from r1 r2
    k_on = r1*r2/k_ph;
    k_off = -(r1+r2)-(k_on+k_ph);
    
    %store info
    sto = [sto; k_on k_off k_ph R_0];
    
    %random initial guesses from output
    r10 = r1*(1 + rjump*(rand-0.5));
    r20 = r2*(1 + rjump*(rand-0.5));
    kph0=k_ph*(1 + rjump*(rand-0.5));
    R00=R_0*(1 + rjump*(rand-0.5));
    
    %decrease annealing temp
    rjump = rjump*0.99;
end

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
% subplot(2,3,6)
plot(t1,R1,'.',t1(end)+dt+t2,R2,'.')
hold on
plot([t1 t1(end)+linspace(0,dt,100) (t2+dt+t1(end))],dub_mol_fun_helper(sol,fitdat),'r')
ylabel('final fit')