% Solve steady state for benchmark and match Kitao
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

%% Initialize
amin = 0.01;
amax = 30;             % clear labor supply
%amax = 800;           % match investment pattern
%amax = 10000;

Na = 1000;
Nk = 1000;
Ne = length(eta);
Nt = length(theta);
N = Na*Ne*Nt;
N_exo = Ne*Nt;

%a = (linspace(0,Na,Na)'/Na).^2*(amax+amin)-amin;
a = linspace(amin,amax,Na)';
aa = repmat(a,[1,N_exo]);
a_s = reshape(aa,N,1);

ee = repmat(eta,[Na,Nt]);
ee_s = reshape(ee,N,1);

tt = repelem(theta,Na,Ne);
tt_s = reshape(tt,N,1);

% iteration parameters
maxIter_vf = 500; 
maxIter_price = 20;       
maxIter_density = 5000;
eps_vf = 1e-7;
eps_price = 1e-4;
eps_density = 1e-10;

%% Iteration on K/L using bisection
%xi = 6.14;      % initial K/L guess, normalize w=1 
%xi = 8.5;
xi = 7.5;        % initial guess
xi_min = 6.5;
xi_max = 8.5;

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
    Iw = wage + r*aa;
    Tw = T(Iw);

    % worker's consumption
    for i = 1:N_exo
        Cw(:,:,i) = bsxfun(@minus,wage(:,i)+assets(:,i)-Tw(:,i),a')./(1+tau_c);
    end
    Cw(Cw<0) = 1e-15;
    Uw = Cw.^(1-sigma)./(1-sigma);

    % solve entrepreneur's profit
    T = @(I) a0.*(I - (I.^(-a1) + a2).^(-1/a1)) + tau_i.*I;
    for i = 1:Na
        for j = 1:N_exo

            ability = theta(floor((j-1)/Ne)+1); 
            if mod(j,Ne) == 0
                units = eta(Ne);
            else
                units = eta(mod(j,Ne));
            end

            kmax = (1+d)*a(i);
            nmax = 30;  % 2000 for unconstrained and amax=800
            kgrid = repmat(linspace(0,kmax,Nk)',1,Nk);  % row vector
            ngrid = repmat(linspace(0,nmax,Nk),Nk,1);   % column vector

            Igrid = max(ability*kgrid.^(alpha*v).*ngrid.^((1-alpha)*v) - delta*kgrid...
                - w*max(ngrid-units,0)...
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
        Ce(:,:,i) = bsxfun(@minus,pi(:,i),a')./(1+tau_c);
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
    inda = discretize(20,a);

    lindw = discretize(Aw,a);
    linde = discretize(Ae,a);
    lindw(Aw==amax) = Na;
    linde(Ae==amax) = Na;
    rindw = lindw + 1;
    rinde = linde + 1;
    rindw(Aw==amax) = Na;
    rinde(Ae==amax) = Na;


    % assign probabilities by lottery
    pleftw = (a(rindw) - Aw)./(a(rindw) - a(lindw));
    pleftw(rindw==lindw) = 1;
    prightw = 1 - pleftw;

    plefte = (a(rinde) - Ae)./(a(rinde) - a(linde));
    plefte(rinde==linde) = 1;
    prighte = 1 - plefte;

    f0 = zeros(Na,N_exo*2);  % initial density
    f0(floor(inda/2),1:N_exo) = 1/N_exo;
    %f0(:,N_exo) = 1/Na;
    %f0(Na/2,N_exo) = 1;
    f1 = zeros(Na,N_exo*2);
    
    % update f1 (32 seconds),distribution not very smooth (a grid)
    for n_density = 1:maxIter_density

        for j = 1:N_exo
            % worker to worker
            f1(lindw(1:inda,j),1:N_exo) = f1(lindw(1:inda,j),1:N_exo) + f0(1:inda,j).*Qs(j,:).*Ow(1:inda,j).*pleftw(1:inda,j);
            f1(rindw(1:inda,j),1:N_exo) = f1(rindw(1:inda,j),1:N_exo) + f0(1:inda,j).*Qs(j,:).*Ow(1:inda,j).*prightw(1:inda,j);

            % worker to entrepreneur
            f1(lindw(1:inda,j),N_exo+1:2*N_exo) = f1(lindw(1:inda,j),N_exo+1:2*N_exo) + f0(1:inda,j).*Qs(j,:).*(1-Ow(1:inda,j)).*pleftw(1:inda,j);
            f1(rindw(1:inda,j),N_exo+1:2*N_exo) = f1(rindw(1:inda,j),N_exo+1:2*N_exo) + f0(1:inda,j).*Qs(j,:).*(1-Ow(1:inda,j)).*prightw(1:inda,j);

            % entrepreneur to worker
            f1(linde(1:inda,j),1:N_exo) = f1(linde(1:inda,j),1:N_exo) + f0(1:inda,j+N_exo).*Qs(j,:).*Oe(1:inda,j).*plefte(1:inda,j);
            f1(rinde(1:inda,j),1:N_exo) = f1(rinde(1:inda,j),1:N_exo) + f0(1:inda,j+N_exo).*Qs(j,:).*Oe(1:inda,j).*prighte(1:inda,j);

            % entrepreneur to entrepreneur
            f1(linde(1:inda,j),N_exo+1:2*N_exo) = f1(linde(1:inda,j),N_exo+1:2*N_exo) + f0(1:inda,j+N_exo).*Qs(j,:).*(1-Oe(1:inda,j)).*plefte(1:inda,j);
            f1(rinde(1:inda,j),N_exo+1:2*N_exo) = f1(rinde(1:inda,j),N_exo+1:2*N_exo) + f0(1:inda,j+N_exo).*Qs(j,:).*(1-Oe(1:inda,j)).*prighte(1:inda,j);
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
Iw = ee*w + aa*r;
Ie = tt.*investment.^(alpha*v).*labor.^((1-alpha)*v) - delta*investment...
    -r*(investment-aa).*(investment<=aa)-(r+iota)*(investment-aa).*(investment>aa)...
    -w*max(labor-ee,0);
income = sum(sum(Iw.*f0(:,1:N_exo) + Ie.*f0(:,N_exo+1:2*N_exo)));
frac_income_e = sum(sum(Ie.*f0(:,N_exo+1:2*N_exo)))/income;

% tax
incometax = sum(sum((T(Iw)-tau_i*Iw).*f0(:,1:N_exo) + (T(Ie)-tau_i*Ie).*f0(:,N_exo+1:2*N_exo)));
ctax_w = (ee*w + (1+r)*aa - T(Iw) - Aw)./(1+tau_c)*tau_c;
ctax_e = (pi - Ae)./(1+tau_c)*tau_c;
tax_total = sum(sum((ctax_w + T(Iw)).*f0(:,1:N_exo)...
    + (ctax_e + T(Ie)).*f0(:,N_exo+1:2*N_exo)));
frac_incometax = incometax / tax_total;
frac_G = tax_total / output;

% wealth distribution 
fa = sum(f0,2);
Fa = cumsum(fa);

p = 0.01:0.01:1;
for pp = 1:100
    [~, pindex(pp)] = min(abs(Fa-p(pp)));
end
perc = a(pindex);

integrand = a.*fa;
top60 = sum(integrand(pindex(40):pindex(end)))/asset_supply*100;
top40 = sum(integrand(pindex(60):pindex(end)))/asset_supply*100;
top20 = sum(integrand(pindex(80):pindex(end)))/asset_supply*100;
top10 = sum(integrand(pindex(90):pindex(end)))/asset_supply*100;
top5 = sum(integrand(pindex(95):pindex(end)))/asset_supply*100;
top1 = sum(integrand(pindex(99):pindex(end)))/asset_supply*100;

S_a = cumsum(a.*fa)/asset_supply;
trapez_a = 0.5*(S_a(1)*fa(1) + sum((S_a(2:Na)+S_a(1:Na-1)).*fa(2:Na)));
gini = 1-2*trapez_a;

% composition of entrepreneurs
t_w = [sum(sum(f0(:,1:5))),sum(sum(f0(:,6:10))),sum(sum(f0(:,11:15))),sum(sum(f0(:,16:20)))];
t_e = [sum(sum(f0(:,21:25))),sum(sum(f0(:,26:30))),sum(sum(f0(:,31:35))),sum(sum(f0(:,36:40)))];
t_all = t_w + t_e;

leverage = max(0,investment./aa-1);
t_asset = f0(:,N_exo+1:N_exo*2).*aa;
t_inv = f0(:,N_exo+1:N_exo*2).*investment;
t_lev = f0(:,N_exo+1:N_exo*2).*leverage;
t_asset = [sum(sum(t_asset(:,1:5))),sum(sum(t_asset(:,6:10))),sum(sum(t_asset(:,11:15))),sum(sum(t_asset(:,16:20)))];
t_inv = [sum(sum(t_inv(:,1:5))),sum(sum(t_inv(:,6:10))),sum(sum(t_inv(:,11:15))),sum(sum(t_inv(:,16:20)))];
t_lev = [sum(sum(t_lev(:,1:5))),sum(sum(t_lev(:,6:10))),sum(sum(t_lev(:,11:15))),sum(sum(t_lev(:,16:20)))];
t_asset = t_asset./t_e;
t_inv = t_inv./t_e;
t_lev = t_lev./t_e;

% asset cdf by worker and entrepreneur
asset_w = sum(sum(f0(:,1:N_exo).*aa));
cdf_w = cumsum(sum(f0(:,1:N_exo),2))/(1-frac_e)*100;
cdf_e = cumsum(sum(f0(:,N_exo+1:N_exo*2),2))/frac_e*100;
pdf_w = sum(f0(:,1:N_exo),2);
pdf_e = sum(f0(:,N_exo+1:N_exo*2),2);

save eqm_kitaoparams_largegrid.mat;

toc;

%% graphs
figure(1)
set(gca,'FontSize',16)
plot(theta,t_all,theta,t_w,theta,t_e)
xlabel('Entrepreneurial ability $\theta$','FontSize',16,'interpreter','latex')
ylabel('Percentage (%)','FontSize',16,'interpreter','latex')

figure(2)
plot(a,investment(:,10),a,investment(:,15),a,investment(:,20),a,a,'--')
xlabel('Assets')
ylabel('Investment')

figure(3)
subplot(1,2,1)
plot(a,pdf_w,a,pdf_e)
xlabel('Assets')
ylabel('Density')

subplot(1,2,2)
plot(a,cdf_w,a,cdf_e)
xlabel('Assets')
ylabel('Percentage (%)')
ylim([0 100])


