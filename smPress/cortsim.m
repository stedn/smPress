% Francois Robin

function [R,Y,NB,NDIS] = cortsim(R0,Y0,b_rate,b_std,dep_rate,ph_rate,syn_rate,simlength,dt)

%Computes sequentially the sequence of values for R, Y, nb (birth) and ndis
%(disappearances) for a time period of simlength (in s), with a time delay
%of dt (eg: .1 for 10 fr/s).

nfr = round(simlength/dt);

R = NaN(1,nfr+1);
Y = NaN(1,nfr+1);
NB = NaN(1,nfr);
NDIS = NaN(1,nfr);

R(1) = R0;
Y(1) = Y0;

dep_rate = dep_rate*dt;
ph_rate = ph_rate*dt;
b_rate = b_rate*dt;
b_std = b_std*dt;
syn_rate = syn_rate*dt;

for frame = 1:nfr
    ndis = poissrnd(R(frame)*(dep_rate+ph_rate));
    ndis = min(R(frame),ndis);
    nd = binornd(ndis, dep_rate/(dep_rate+ph_rate));
    np = ndis-nd;
    nb = poissrnd(b_rate*Y(frame));
    NB(frame) = nb;
    NDIS(frame) = ndis;
    Y(frame+1) = Y(frame)-nb+nd+poissrnd(syn_rate);
    R(frame+1) = R(frame)+nb-nd-np;
end

end

