function dy = nmol_ode(t,y,kon, koff, kph, ks)
    dy = [kon*y(2) - (koff+kph)*y(1); -kon*y(2) + koff*y(1) + ks];
end