% Will McFadden (wmcfadden) and Francois Robin
% a test to see how varying the number of runs and duty ratio of the data
% has an impact on measurement accuracy for data simulated with cortsim.

%% Setting the parameters

exp_time = 200;     %approximate max duration of maintenance phase 
nsim = 10;         %number of simulations to run

itn = 4;            %number of iterations of the on-off cycle

% duty = .5;          %fraction of the time with laser "on"
stordif = [];
storerr = [];
cc = jet(8*9);
cci = 1;
for ph_rate = linspace(0.05,0.5,10)
    for duty = 0.1:0.1:0.9
        dt = 0.1;
        cycle_time = exp_time/(itn-1+duty);
        cycle_steps = round(cycle_time/dt);   %number of steps in a cycle
        on_steps = round(duty*cycle_steps);
        off_steps = cycle_steps-on_steps;
        tmax = on_steps*dt;
%         on_steps = tmax/dt;
        twait = off_steps*dt;
        actual_duration = (twait+tmax)*(itn-1)+tmax;
        total_steps = actual_duration/dt+1;
        unused_time = exp_time - actual_duration;

        b_rate = 0.01;	% probability of birth
        syn_rate = 1;       % number of molecules synthesized per second
        dep_rate = .1;   % depolymerization rate (only one? - is there a reliable way to say if a better mix would be closer as well as inexpensive in hypotheses?)
%         ph_rate= .1;		% photobleaching rate
        R0_ = 250;       % Number of molecules at the cortex at t=0
        b_std = b_rate/5;
        
        Y0_ = round(R0_*dep_rate/b_rate);         % number of molecules in the cytoplasm
        N0_ = R0_+Y0_;		% number of molecules on the cortex at t0

        Rg = NaN(on_steps+1,itn,nsim);
        Yg = NaN(on_steps+1,itn,nsim);
        NBg = NaN(on_steps,itn,nsim);

        for i = 1:nsim

            %% Defining output variables

            R0 = R0_;
            N0 = N0_;
            Y0 = Y0_;

            R = NaN(on_steps+1,itn);
            Y = NaN(on_steps+1,itn);
            NB = NaN(on_steps,itn);
            NDIS = NaN(on_steps,itn);
            T = [0:dt:tmax]'*ones(1,itn);
            DT = twait*ones(1,itn-1);


            %% First cycle, not looped

                % Randomization of R0 and Y0, run for 30 cycles (convergence is very rapid, 5 cycles should be enough, but just in case), separated by 1s
            [R0est,Y0est,NBest,NDISest] = cortsim(R0,Y0,b_rate,b_std,dep_rate,0,0,3/dep_rate,dt);

            R(1,1) = R0est(1,end);		% number of molecules in the cytoplasm
            Y(1,1) = Y0est(1,end);		% number of molecules on the cortex at t0

                % First bleaching phase

            Rint = [];
            Yint = [];
            NBint = [];

            [Rint,Yint,NBint] = cortsim(R(1,1),Y(1,1),b_rate,b_std,dep_rate,ph_rate,syn_rate,tmax,dt);
            NB(:,1) = NBint';
            Y(:,1) = Yint';
            R(:,1) = Rint';

            %% Subsequent cycles, looped

            for j = 2:itn
                Rint = [];
                Yint = [];
                NBint = [];
                R0 = R(end,j-1);
                Y0 = Y(end,j-1);
                [Rint,Yint,NBint] = cortsim(R0,Y0,b_rate,b_std,dep_rate,0,syn_rate,twait,dt);

                R0 = Rint(1,end);
                Y0 = Yint(1,end);

                Rint = [];
                Yint = [];
                NBint = [];
                [Rint,Yint,NBint] = cortsim(R0,Y0,b_rate,b_std,dep_rate,ph_rate,syn_rate,tmax,dt);
                R(:,j) = Rint';
                Y(:,j) = Yint';    
                NB(:,j) = NBint';
            end
            Rg(:,:,i) = R;
            Yg(:,:,i) = Y;
            NBg(:,:,i) = NB;
        end

        R = mean(Rg,3);
        NB = mean(NBg,3);

        fitdat = {T DT};
        fitval = [R(:)' NB(:)'];

        kont = b_rate;
        kofft = dep_rate;
        kpht = ph_rate;
        kst = syn_rate;
        R0t = R(1,1);
        r1t = -(kont+kofft+kpht)/2 - sqrt(((kont+kofft+kpht))^2-4*kont*kpht)/2;
        r2t = -(kont+kofft+kpht)/2 + sqrt(((kont+kofft+kpht))^2-4*kont*kpht)/2;


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
            
            weights = 1./fitval;
            weight_fun = @(q,dat) weights .* n_mol_fun(q,dat);
            [sol,MSE,residual,exitflag,output,lambda,J] = lsqcurvefit(weight_fun,param,fitdat,weights.*fitval,lb,up);

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
        MSE = MSE/length(fitval);
        covB = inv(J'*J)*MSE;
        j_on = [r2/k_ph r1/k_ph -r1*r2/k_ph];
        j_off = [1+r2/k_ph 1+r1/k_ph 1 - r1*r2/k_ph];
        subcov = covB(1:3,1:3);
        var_kon = j_on*subcov*j_on';
        var_koff = j_off*subcov*j_off';
        var_kph = covB(3,3);
        var_ks = covB(4,4);
        var_R0 = covB(5,5);

        t = T;
        dt = DT;
        fsol = n_mol_fun_full(sol,fitdat);
        R_sol = fsol(1:end/2);
        birth_sol = fsol(end/2+1:end);
        figure;
        subplot(2,1,1);

        offsett = 0;

        for i=1:length(dt)
            plot(offsett + t(:,i),R(:,i),'.')
            offsett = offsett + t(end,i) + dt(i);
            hold on
        end
        plot(offsett + t(:,end),R(:,end),'.')
        plot((1:length(R_sol))*(t(2,1)-t(1,1)), R_sol,'r');  
        ylabel('final fit to R')

        offsett = 0;


        subplot(2,1,2)

        for i=1:length(dt)
            plot(offsett + t(2:end,i),NB(:,i),'.')
            offsett = offsett + t(end,i) + dt(i);
            hold on
        end
        plot(offsett + t(2:end,end),NB(:,end),'.')
        plot((1:length(birth_sol))*(t(2,1)-t(1,1)), birth_sol,'r');  

        ylabel('final fit to births')
        
        finalanswer = full([[k_on; k_off; k_ph; k_s; R_0] sqrt([var_kon;var_koff;var_kph;var_ks;var_R0])])
        francsanswer = [kont; kofft; kpht; kst; R0t]
        stordif = [stordif [k_on; k_off; k_ph; k_s; R_0]-[kont; kofft; kpht; kst; R0t]];
        storerr = [storerr sqrt([var_kon;var_koff;var_kph;var_ks;var_R0])];
        cci = cci +1;
    end 
end
toshow = abs(stordif(1,:));
[X,Y] = meshgrid(linspace(0.05,0.5,10),0.1:0.1:0.9);
Z = reshape(toshow,size(X,1),size(X,2));
surf(X,Y,Z);
