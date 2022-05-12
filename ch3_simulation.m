% Simulate lottery sector
clear all; close all; clc;
% load('eqm_lottery_largegrid.mat');
% load('eqm_lottery_small_largegrid.mat');
 load('eqm_lottery_large_largegrid.mat');
tic;

%% Initialize 
I = 400000;                     % household sample size
T_use = 100;                    % useful time periods, use t=199
T_burn = 100;                   % burning periods
T = T_use+T_burn;

wealth = zeros(I,T);            % start at steady state capital
consumption = zeros(I,T);
occupation = zeros(I,T);        % occupation, =1 if worker, half
capital_choice = zeros(I,T);

%% Generate shocks
% shocks for eta
shock_eta = zeros(I,T);
shock_eta_val = zeros(I,T);

shock_eta_val(:,1) = randsample(eta,I,true,st_eta);
shock_eta(:,1) = discretize(shock_eta_val(:,1),eta);
shock_eta(shock_eta_val(:,1)==eta(Ne),1) = Ne;

% shocks for theta
shock_theta = zeros(I,T);
shock_theta_val = zeros(I,T);

shock_theta_val(:,1) = randsample(theta,I,true,st_theta);
shock_theta(:,1) = discretize(shock_theta_val(:,1),theta);
shock_theta(shock_theta_val(:,1)==theta(Nt),1) = Nt;

% shocks for lottery prizes
prize_amount = randsample(vprize,I*T,true,vprob);
prize_ind = discretize(prize_amount,vprize);
prize_ind(prize_amount==vprize(Ny)) = Ny;
prize_amount = reshape(prize_amount,I,T);
prize_ind = reshape(prize_ind,I,T);



%% Simulate and update policy functions
% known: c_e, c_w, Ae, Aw, Oe, Ow
% move forward: wealth, consumption, occupation, investment
wealth(:,1) = asset_supply;
%occupation(:,1) = randsample([0,1],I,true,[0.5,0.5]);
occupation(:,1) = 1; % prevent theta=0 entrepreneurs

for t = 2:T
    if mod(t,10) == 0
        disp('updating position at')
        disp(t)
    end
    
    for i = 1:I
        % update eta and theta
        shock_eta(i,t) = randsample(Ne,1,true,PI_eta(shock_eta(i,t-1),:));
        shock_theta(i,t) = randsample(Nt,1,true,PI_theta(shock_theta(i,t-1),:));

        shock_eta_val(i,t) = eta(shock_eta(i,t));
        shock_theta_val(i,t) = theta(shock_theta(i,t));

        index = (prize_ind(i,t)-1)*Ne*Nt + (shock_theta(i,t)-1)*Ne + shock_eta(i,t);

        if occupation(i,t-1) == 1 % worker
            % update wealth, consumption, occupation
            wealth(i,t) = max(0,interp1(a,Aw(:,index),wealth(i,t-1),'linear','extrap'));
            consumption(i,t) = max(0,interp1(a,c_w(:,index),wealth(i,t-1),'linear','extrap'));
            occupation(i,t) = interp1(a,Ow(:,index),wealth(i,t-1),'linear','extrap');
            occupation(i,t) = 0*(occupation(i,t)<0.5) + 1*(occupation(i,t)>=0.5);

        else
            % update wealth, consumption, occupation, capital choice
            wealth(i,t) = max(0,interp1(a,Ae(:,index),wealth(i,t-1),'linear','extrap'));
            consumption(i,t) = max(0,interp1(a,c_e(:,index),wealth(i,t-1),'linear','extrap'));
            occupation(i,t) = interp1(a,Oe(:,index),wealth(i,t-1),'linear','extrap');
            occupation(i,t) = 0*(occupation(i,t)<0.5) + 1*(occupation(i,t)>=0.5);
            capital_choice(i,t) = max(0,interp1(a,investment(:,index),wealth(i,t-1),'linear','extrap'));

        end   

    end

end


%save simulation_fourprize.mat;

reg_matrix = [occupation(:,T-1),occupation(:,T),wealth(:,T-1),wealth(:,T),consumption(:,T-1),consumption(:,T),...
    capital_choice(:,T-1),capital_choice(:,T),prize_amount(:,T-1),prize_amount(:,T),shock_theta_val(:,T-1),shock_theta_val(:,T)];
writematrix(reg_matrix,'reg_matrix_large_newgrid.csv');

toc;