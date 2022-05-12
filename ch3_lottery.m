% Solve steady state for lottery sector - 13845s
clear all; close all; clc;

%% Parameter
% Preference
sigma = 2.0;                            % literature
beta = 0.9575;                          % discount factor
%beta = 0.9428;

% Production technology
v = 0.88;                               % non-corporate production parameter
alpha = 0.36;                           % capital share in the corporate sector
delta = 0.06;                           % depreciation rate of capital
%A = 0.8129;                            % normalize w = 1
A = 1;

% Intermediary sector
d = 0.50;                               % maximum leverage ratio
iota = 0.05;                            % loan premium

% lottery
phi = 0.0286;                           % lottery ticket price
gamble_perc = 1.32/100;

size_ratio = [0,1,3,6];
vprob = [0.9950, 0.0047, 0.00025, 0.00005];
prize_small = phi/sum(size_ratio.*vprob);
vprize = prize_small*size_ratio;

% Tax
tau_c = 0.0567;                         % consumption tax
tau_i = 0.0316;                         % income tax rate
a0 = 0.258;
a1 = 0.768;
a2 = 0.438;
T = @(I) a0.*(I - (I.^(-a1) + a2).^(-1/a1)) + tau_i.*I;

% Productivity and entrepreneurial ability
eta = [0.646, 0.798, 0.966, 1.169, 1.444];
PI_eta = [ 0.731, 0.253, 0.016, 0.000, 0.000;
           0.192, 0.555, 0.236, 0.017, 0.000;
           0.011, 0.222, 0.534, 0.222, 0.011; % modeified 0.534 to add up to 1
           0.000, 0.017, 0.236, 0.555, 0.192;
           0.000, 0.000, 0.016, 0.253, 0.731];
[st_eta x] = eigs(PI_eta',1);
st_eta = st_eta./sum(st_eta);           % stationary distribution for eta

%tbar = 1.55;      % scaling factor (match investment policy)
tbar = 1;
theta = tbar*[0.000, 0.706, 1.470, 2.234];
PI_theta = [ 0.780, 0.220, 0.000, 0.000;
             0.430, 0.420, 0.150, 0.000;
             0.000, 0.430, 0.420, 0.150;
             0.000, 0.000, 0.220, 0.780];
[st_theta x] = eigs(PI_theta',1);
st_theta = st_theta./sum(st_theta);     % stationary distribution for eta

Qs = kron(PI_theta,PI_eta);
Qs = kron(repmat(vprob,length(vprize),1),Qs);

%% Initialize
amin = 0.01;
amax = 20;             % clear labor supply
%amax = 800;           % match investment pattern
%amax = 10000;

Na = 1000;
Nk = 500;
Ne = length(eta);
Nt = length(theta);
Nl = 1000;
Ny = length(vprize);
N = Na*Ne*Nt*Ny;
N_exo = Ne*Nt*Ny;

%a = (linspace(0,Na,Na)'/Na).^2*(amax+amin)-amin;
a = linspace(amin,amax,Na)';
afiner = linspace(amin,amax,Nl)';
inda = discretize(a,afiner);
inda(a==amax) = Nl;
aa = repmat(a,[1,N_exo]);
a_s = reshape(aa,N,1);

ee = repmat(eta,[Na,Nt*Ny]);
ee_s = reshape(ee,N,1);

tt = repelem(theta,Na,Ne);
tt = repmat(tt,1,Ny);
tt_s = reshape(tt,N,1);

yy = repelem(vprize,Na,Ne*Nt);
yy_s = reshape(yy,N,1);

% iteration parameters
maxIter_vf = 500; 
maxIter_price = 15;  
maxIter_density = 5000;
eps_vf = 1e-7;
eps_price = 1e-4;
eps_density = 1e-10;

%% Iteration on K/L using bisection
%xi = 6.14;      % initial K/L guess, normalize w=1 
%xi = 8.5;
xi = 7.25;        % initial guess
xi_min = 6.5;
xi_max = 8.0;

tic;
for n_price = 1:maxIter_price
    % start loop for K/L here
    r = A*alpha*xi^(alpha-1) - delta;
    w = A*(1-alpha)*xi^alpha;

    % value function iteration (E,W) - solve consumption, savings, and
    % occupation
    v0e = zeros(Na,N_exo);
    v0w = zeros(Na,N_exo);
    v1w = zeros(Na,N_exo);
    v1e = zeros(Na,N_exo);
    v1ww = zeros(Na,N_exo);
    v1we = zeros(Na,N_exo);
    v1ew = zeros(Na,N_exo);
    v1ee = zeros(Na,N_exo);
    Ce = zeros(Na,Na,N_exo);
    Cw = zeros(Na,Na,N_exo);
    Ae = zeros(Na,N_exo);
    Aw = zeros(Na,N_exo);
    Oe = zeros(Na,N_exo);
    Ow = zeros(Na,N_exo);
    pi = zeros(Na,N_exo);
    investment = zeros(Na,N_exo);
    labor = zeros(Na,N_exo);

    % update consumption using budget constraint
    wage = ee*w;
    assets = (1+r)*aa;
    Iw = wage + r*aa + yy - phi;
    Tw = T(Iw);

    % worker's consumption
    for i = 1:N_exo
        Cw(:,:,i) = bsxfun(@minus,wage(:,i)+assets(:,i)+yy(:,i)-phi-Tw(:,i),a')./(1+tau_c);
    end
    Cw(Cw<0) = 1e-15;
    Uw = Cw.^(1-sigma)./(1-sigma);

    % solve entrepreneur's profit
    T = @(I) a0.*(I - (I.^(-a1) + a2).^(-1/a1)) + tau_i.*I;
    for i = 1:Na
        for j = 1:N_exo
            
            ability = tt(1,j);
            units = ee(1,j);
            prize = yy(1,j);
            
            kmax = (1+d)*(a(i)+prize);
            nmax = 30;  % 2000 for unconstrained and amax=800
            kgrid = repmat(linspace(0,kmax,Nk)',1,Nk);  % row vector
            ngrid = repmat(linspace(0,nmax,Nk),Nk,1);   % column vector

            Igrid = max(ability*kgrid.^(alpha*v).*ngrid.^((1-alpha)*v) - delta*kgrid...
                - w*max(ngrid-units,0) - phi + prize...
                - r*(kgrid-a(i)).*(kgrid<a(i)) - (r+iota)*(kgrid-a(i)).*(kgrid>a(i)),0); 

            Tgrid = T(Igrid);

            pigrid = ability*kgrid.^(alpha*v).*ngrid.^((1-alpha)*v) + (1-delta)*kgrid...
                - (1+r)*(kgrid-a(i)).*(kgrid<a(i)) - (1+r+iota)*(kgrid-a(i)).*(kgrid>a(i))...
                - w*max(ngrid-units,0) - Tgrid;

            maxpi = max(max(pigrid));
            [k,n] = find(pigrid == maxpi);
            pi(i,j) = maxpi;
            investment(i,j) = kgrid(k(1));
            labor(i,j) = ngrid(1,n(1));

        end
    end

    % entrepreneur's consumption
    for i = 1:N_exo
        Ce(:,:,i) = bsxfun(@minus,pi(:,i)+yy(:,i)-phi,a')./(1+tau_c);
    end
    Ce(Ce<0) = 1e-15;
    Ue = Ce.^(1-sigma)./(1-sigma);

    % update both value functions, loop starts here
    for n_vf = 1:maxIter_vf
        EVw = v0w*Qs';
        EVe = v0e*Qs';

        for j = 1:N_exo
            [v1ww(:,j),indw] = max(Uw(:,:,j) + beta*repmat(EVw(:,j)',Na,1),[],2);
            [v1we(:,j),inde] = max(Uw(:,:,j) + beta*repmat(EVe(:,j)',Na,1),[],2);

            v1w(:,j) = max(v1ww(:,j),v1we(:,j));
            Ow(:,j) = 1.0.*(v1ww(:,j) >= v1we(:,j));  % 1 if worker
            Aw(:,j) = a(indw).*(v1ww(:,j) >= v1we(:,j)) + a(inde).*(v1ww(:,j) < v1we(:,j));

            [v1ew(:,j),indw] = max(Ue(:,:,j) + beta*repmat(EVw(:,j)',Na,1),[],2);
            [v1ee(:,j),inde] = max(Ue(:,:,j) + beta*repmat(EVe(:,j)',Na,1),[],2);

            v1e(:,j) = max(v1ew(:,j),v1ee(:,j));
            Oe(:,j) = 1.0.*(v1ew(:,j) >= v1ee(:,j));
            Ae(:,j) = a(indw).*(v1ew(:,j) >= v1ee(:,j)) + a(inde).*(v1ew(:,j) < v1ee(:,j));

        end

        error = max(max(max(abs(v1w-v0w))), max(max(abs(v1e-v0e))));

        if error < eps_vf
            break
        else
            v0w = v1w;
            v0e = v1e;
        end
    end

    % save vfibenchmark_smallkgrid.mat;

    % find invariant distribution
    % find index of policy a in finer grid of a
    lindw = discretize(Aw,afiner);
    linde = discretize(Ae,afiner);
    lindw(Aw==amax) = Nl;
    linde(Ae==amax) = Nl;
    rindw = lindw + 1;
    rinde = linde + 1;
    rindw(Aw==amax) = Nl;
    rinde(Ae==amax) = Nl;


    % assign probabilities by lottery
    pleftw = (afiner(rindw) - Aw)./(afiner(rindw) - afiner(lindw));
    pleftw(rindw==lindw) = 1;
    prightw = 1 - pleftw;

    plefte = (afiner(rinde) - Ae)./(afiner(rinde) - afiner(linde));
    plefte(rinde==linde) = 1;
    prighte = 1 - plefte;

    f0 = zeros(Nl,N_exo*2);  % initial density
    %f0(Nl,1:N_exo) = 1/N_exo;
    f0(:,N_exo) = 1/Na;
    %f0(Na/2,N_exo) = 1;
    f1 = zeros(Nl,N_exo*2);

    % update f1 (32 seconds),distribution not very smooth (a grid)
    for n_density = 1:maxIter_density

        for j = 1:N_exo
            % worker to worker
            f1(lindw(:,j),1:N_exo) = f1(lindw(:,j),1:N_exo) + f0(inda,j).*Qs(j,:).*Ow(:,j).*pleftw(:,j);
            f1(rindw(:,j),1:N_exo) = f1(rindw(:,j),1:N_exo) + f0(inda,j).*Qs(j,:).*Ow(:,j).*prightw(:,j);

            % worker to entrepreneur
            f1(lindw(:,j),N_exo+1:2*N_exo) = f1(lindw(:,j),N_exo+1:2*N_exo) + f0(inda,j).*Qs(j,:).*(1-Ow(:,j)).*pleftw(:,j);
            f1(rindw(:,j),N_exo+1:2*N_exo) = f1(rindw(:,j),N_exo+1:2*N_exo) + f0(inda,j).*Qs(j,:).*(1-Ow(:,j)).*prightw(:,j);

            % entrepreneur to worker
            f1(linde(:,j),1:N_exo) = f1(linde(:,j),1:N_exo) + f0(inda,j+N_exo).*Qs(j,:).*Oe(:,j).*plefte(:,j);
            f1(rinde(:,j),1:N_exo) = f1(rinde(:,j),1:N_exo) + f0(inda,j+N_exo).*Qs(j,:).*Oe(:,j).*prighte(:,j);

            % entrepreneur to entrepreneur
            f1(linde(:,j),N_exo+1:2*N_exo) = f1(linde(:,j),N_exo+1:2*N_exo) + f0(inda,j+N_exo).*Qs(j,:).*(1-Oe(:,j)).*plefte(:,j);
            f1(rinde(:,j),N_exo+1:2*N_exo) = f1(rinde(:,j),N_exo+1:2*N_exo) + f0(inda,j+N_exo).*Qs(j,:).*(1-Oe(:,j)).*prighte(:,j);
        end

        f1 = f1./(sum(sum(f1)));    

        if max(max(abs(f1-f0))) < eps_density
            break
        end

        f0 = f1;
    end

    %save dist_kitaoparams_smalla.mat

    % update capital-labor ratio
    asset_supply = sum(sum(f0.*repmat(aa,1,2)));
    investment_e = sum(sum(f0(:,N_exo+1:2*N_exo).*investment));
    investment_c = asset_supply - investment_e;

    labor_supply = sum(sum(f0.*repmat(ee,1,2)));
    labor_e = sum(sum(f0(:,N_exo+1:2*N_exo).*labor));
    labor_c = labor_supply - labor_e;

    xi_new = investment_c/labor_c;

    if xi_new - xi > eps_price
        xi_min = xi;
        xi = 0.5*(xi_min + xi_max);
    elseif xi_new - xi < -eps_price
        xi_max = xi;
        xi = 0.5*(xi_min + xi_max);
    else
        break
    end
    
end

% checks
total_prize = sum(sum((f0(:,1:80)+f0(:,81:160)).*yy));

%% aggregate quantities
frac_e = sum(sum(f0(:,N_exo+1:2*N_exo)));
exit_rate = sum(sum(f0(:,N_exo+1:2*N_exo).*Oe))/sum(sum(f0(:,N_exo+1:2*N_exo)));
%sum(sum(f0(:,N_exo+1:2*N_exo).*Oe))/sum(sum(f0(:,N_exo+1:2*N_exo)+f0(:,1:N_exo).*(1-Ow)));
frac_capital_e = investment_e / asset_supply;
asset_e = sum(sum(f0(:,N_exo+1:2*N_exo).*aa));
frac_asset_e = asset_e / asset_supply;

output_c = A*investment_c^alpha*labor_c^(1-alpha);
output_e = sum(sum(tt.*investment.^(alpha*v).*labor.^((1-alpha)*v).*f0(:,N_exo+1:2*N_exo)));
output = output_c + output_e;
KYratio = asset_supply/output;

% income
Iw = ee*w + aa*r + yy - phi;
Ie = tt.*investment.^(alpha*v).*labor.^((1-alpha)*v) - delta*investment...
    -r*(investment-aa).*(investment<=aa)-(r+iota)*(investment-aa).*(investment>aa)...
    -w*max(labor-ee,0)+yy-phi;
income = sum(sum(Iw.*f0(:,1:N_exo) + Ie.*f0(:,N_exo+1:2*N_exo)));
frac_income_e = sum(sum(Ie.*f0(:,N_exo+1:2*N_exo)))/income;

% consumption
c_w = (ee*w + (1+r)*aa + yy - phi - T(Iw) - Aw)./(1+tau_c);
c_e = (pi + yy - phi - Ae)./(1+tau_c);

% tax
incometax = sum(sum((T(Iw)-tau_i*Iw).*f0(:,1:N_exo) + (T(Ie)-tau_i*Ie).*f0(:,N_exo+1:2*N_exo)));
ctax_w = c_w*tau_c;
ctax_e = c_e*tau_c;
tax_total = sum(sum((ctax_w + T(Iw)).*f0(:,1:N_exo)...
    + (ctax_e + T(Ie)).*f0(:,N_exo+1:2*N_exo)));
frac_incometax = incometax / tax_total;
frac_G = tax_total / output;

% asset cdf by worker and entrepreneur
asset_w = sum(sum(f0(:,1:N_exo).*aa));
cdf_w = cumsum(sum(f0(:,1:N_exo),2))/(1-frac_e)*100;
cdf_e = cumsum(sum(f0(:,N_exo+1:N_exo*2),2))/frac_e*100;
pdf_w = sum(f0(:,1:N_exo),2);
pdf_e = sum(f0(:,N_exo+1:N_exo*2),2);

save eqm_lottery.mat;

toc;