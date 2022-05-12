% General Equilibrium VC effort
% 40 seconds GE loop
clc; clear all; close all; 

%% Parameters
sigma = 5; % vc and founder elasticity
theta = 0.1; % management elasticity is 0.1 (increase to 0.5 also works)
beta = 0.6; % labor elasticity
chi = 1-beta-theta; 
alpha = 0.7; % bargaining power, >0.5 gives cutoff
nu = 0.57; % project success rate
se = 0.96; % firm survival rate 
sv = 0.8; % vc survival rate

gamma = 1; % relative importance of VC
I = 30; % fixed cost, affect z cutoff
ke = 10; % entry cost constant entrepreneur
kv = 10; % entry cost constant vc
etae = 2; % entry cost elasticity entrepreneur
etav = 2; % entry cost elasticity vc
L = 50; % total labor, exogenous

%% Initialize
zmin = 10;
zmax = 100;
cmin = 0.1;
cmax = 3; 
N_z = 100;
N_c = 20;
z = linspace(zmin,zmax,N_z);
c = linspace(cmin,cmax,N_c);
zvec = repmat(z',1,N_c);
cvec = repmat(c,N_z,1);
maxIter = 100;
eps = 1e-3;
relax_m = 0.2;    
relax_h = 0.2;
options = optimset('Display','iter','FunValCheck','on','MaxFunEvals',100000,'MaxIter',1000,'TolX',1e-5,'TolFun',1e-10);

% guess interest rate 
r = 0.05;
%r = 0.0502; % I = 20
delta = 1/(1+r);

% guess eqm mass
M = 10;
H = 1;
    
for hloop = 1:maxIter
    % guess wage and vc cost
    w = 1;
    wmin = 0.1;
    wmax = 5;
    
    for wloop = 1:maxIter
    
        omega = 1;    
        omin = 0.1;
        omax = 10;
        
        for oloop = 1:maxIter
            % Given price, solve for value function and funding status
            % solo
            f_solo = (beta/w)^(beta/chi)*zvec.*(delta*nu*theta/(1-delta*se)./cvec).^((theta+chi)/chi);
            l_solo = (beta/w)^(1/(1-beta)).*zvec.^(chi/(1-beta)).*f_solo.^(theta/(1-beta));
            V_solo = chi*(beta/w)^(beta/chi)*(delta*nu/(1-delta*se))^((theta+chi)/chi)...
                .*(theta./cvec).^(theta/chi).*zvec - I;
            E = max(V_solo,0);
            funded = (V_solo > 0)*1.0;
            
            % vc 
            h_v = (beta/w)^(beta/chi)*gamma^(theta*sigma/chi/(sigma-1))*(alpha*delta*nu*theta/omega/(1-delta*se))^((theta+chi)/chi)...
                .*zvec.*(1+gamma^(-sigma).*(omega*(1-alpha)/alpha./cvec).^(sigma-1)).^(theta/chi/(sigma-1)-1);
            f_v = (omega./cvec*(1-alpha)/alpha/gamma).^sigma.*h_v;
            F_v = (gamma*h_v.^(1-1/sigma)+f_v.^(1-1/sigma)).^(sigma/(sigma-1));
            l_v = (beta/w)^(1/(1-beta)).*zvec.^(chi/(1-beta)).*F_v.^(theta/(1-beta));
            S_e = alpha.*(delta*nu/(1-delta*se)*(1-beta)*(beta/w)^(beta/(1-beta))...
                .*zvec.^(chi/(1-beta)).*F_v.^(theta/(1-beta)) - I - E);
            S_v = (1-alpha).*(delta*nu/(1-delta*se)*(1-beta)*(beta/w)^(beta/(1-beta))...
                .*zvec.^(chi/(1-beta)).*F_v.^(theta/(1-beta)) - I - E);
            vcfunded = (S_e > 0).*1.0;
            
            solo = funded.*(1-vcfunded);
            
            % update total supply of H 
            if sum(sum(vcfunded)) > 0
                h_d = M*(1-se)/nu*sum(sum(h_v.*vcfunded))/N_z/N_c;
            else
                h_d = 0;
            end
        
            % update guess for omega
            if h_d - H > eps
                omin = omega;
                omega = (omin + omax)/2;
            elseif h_d - H < -eps
                omax = omega;
                omega = (omin + omax)/2;
            else
                break
            end
        end
    
        l_d = M*sum(sum(l_solo.*solo+l_v.*vcfunded))/N_z/N_c;

        if l_d - L > eps
            wmin = w;
            w = (wmin + wmax)/2;
        elseif l_d - L < -eps
            wmax = w;
            w = (wmin + wmax)/2;
        else
            break
        end
    end
    
    % update expected value and mass for H
    ev_v = sum(sum(S_v.*vcfunded))/N_c/N_z;
    implied_h = (ev_v/kv)^(1/etav)/(1-sv);
    err_h = max(abs(H-implied_h));
    H = relax_h*implied_h + (1-relax_h)*H;
    
    % update expected value and mass for M
    ev_e = sum(sum(S_e.*vcfunded+V_solo.*funded))/N_c/N_z;
    implied_m = (ev_e/ke)^(1/etae)/(1-se)*nu;
    err_m = max(abs(M-implied_m));
    M = relax_m*implied_m + (1-relax_m)*M;

    disp(M)
    disp(H)
       
    if err_m && err_h < eps
        break
    end

end

% labor productivity (should equal to w/beta)
y_solo = zvec.^chi.*f_solo.^theta.*l_solo.^beta.*solo; 
y_v = zvec.^chi.*F_v.^theta.*l_v.^beta.*vcfunded;
Y = M*sum(sum(y_solo+y_v))/N_z/N_c; % multiply by M here
tfp = Y/L;

%% Plots
load("supplyv4.mat");

figure(1)
plot(z,h_v(:,N_c/2).*vcfunded(:,N_c/2))
xlabel('Total Factor Productivity $z$','interpreter','latex')
ylabel('Optimal Effort $h$','interpreter','latex')
set(gca,'XTick',[], 'YTick', [])

plot(c,h_v(N_z,:).*vcfunded(N_z,:))
xlabel('Founder''s Effort Cost $c$','interpreter','latex')
ylabel('Optimal Effort $h$','interpreter','latex')
set(gca,'XTick',[], 'YTick', [])

figure(2)
pcolor(z,c,solo')
xlabel('Idea Productivity','interpreter','latex')
ylabel('Founder Productivity','interpreter','latex')

figure(3)
pcolor(z,c,vcfunded')
xlabel('Idea Productivity','interpreter','latex')
ylabel('Founder Productivity','interpreter','latex')

figure(4)
plot(z,h_v(:,N_c/2).*vcfunded(:,N_c/2),z,h_v_I(:,N_c/2).*vcfunded_I(:,N_c/2))
xlabel('Idea Productivity $z$','interpreter','latex')
ylabel('Optimal Effort $h$','interpreter','latex')
title('Optimal Effort Levels','interpreter','latex')
legend('Benchmark','Lower $I$','Location','northwest','interpreter','latex')
set(gca,'XTick',[], 'YTick', [])

figure(5)
subplot(1,2,1)
h=pcolor(z,c,vcfunded')
xlabel('Idea Productivity $z$','interpreter','latex')
ylabel('Founder Cost $c$','interpreter','latex')
title('Benchmark','interpreter','latex')
set(gca,'XTick',[], 'YTick', [])
set(h, 'EdgeColor', 'none')

subplot(1,2,2)
h=pcolor(z,c,vcfunded_I')
xlabel('Idea Productivity $z$','interpreter','latex')
ylabel('Founder Cost $c$','interpreter','latex')
title('Lower $I$','interpreter','latex')
set(gca,'XTick',[], 'YTick', [])
set(h, 'EdgeColor', 'none')

figure(6)
plot(z,h_v(:,N_c/2).*vcfunded(:,N_c/2),z,h_v_r(:,N_c/2).*vcfunded_r(:,N_c/2))
xlabel('Idea Productivity','interpreter','latex')
ylabel('Optimal Effort','interpreter','latex')
legend('Benchmark','Lower $r$','interpreter','latex')

figure(7)
subplot(1,2,1)
pcolor(z,c,vcfunded')
xlabel('Idea Productivity','interpreter','latex')
ylabel('Founder Productivity','interpreter','latex')
title('Benchmark','interpreter','latex')

subplot(1,2,2)
pcolor(z,c,vcfunded_r')
xlabel('Idea Productivity','interpreter','latex')
ylabel('Founder Productivity','interpreter','latex')
title('Lower $r$','interpreter','latex')

figure(8)
subplot(1,2,1)
plot(z,h_v(:,N_c/2).*vcfunded(:,N_c/2),z,h_v_I(:,N_c/2).*vcfunded_I(:,N_c/2))
xlabel('Idea Productivity $z$','interpreter','latex')
ylabel('$h$','interpreter','latex')
title('Optimal Effort Levels','interpreter','latex')
legend('Benchmark','Lower $I$','Location','northwest','interpreter','latex')
set(gca,'XTick',[], 'YTick', [])

subplot(1,2,2)
plot(z,S_e(:,N_c/2).*vcfunded(:,N_c/2),z,S_e_I(:,N_c/2).*vcfunded_I(:,N_c/2))
xlabel('Idea Productivity $z$','interpreter','latex')
ylabel('$S_e$','interpreter','latex')
title('Entrepreneur''s Surplus','interpreter','latex')
legend('Benchmark','Lower $I$','Location','northwest','interpreter','latex')
set(gca,'XTick',[], 'YTick', [])

figure(9)
subplot(1,2,1)
plot(z,h_v(:,N_c/2).*vcfunded(:,N_c/2),z,h_v_r(:,N_c/2).*vcfunded_r(:,N_c/2))
xlabel('Idea Productivity $z$','interpreter','latex')
ylabel('$h$','interpreter','latex')
title('Optimal Effort Levels','interpreter','latex')
legend('Benchmark','Lower $r$','Location','northwest','interpreter','latex')
set(gca,'XTick',[], 'YTick', [])

subplot(1,2,2)
plot(z,S_e(:,N_c/2).*vcfunded(:,N_c/2),z,S_e_r(:,N_c/2).*vcfunded_r(:,N_c/2))
xlabel('Idea Productivity $z$','interpreter','latex')
ylabel('$S_e$','interpreter','latex')
title('Entrepreneur''s Surplus','interpreter','latex')
legend('Benchmark','Lower $r$','Location','northwest','interpreter','latex')
set(gca,'XTick',[], 'YTick', [])

figure(10)
subplot(1,2,1)
plot(c,h_v(67,:).*vcfunded(67,:),c,h_v_c(67,:).*vcfunded_c(67,:))
xlabel('Founder''s cost $c$','interpreter','latex')
ylabel('$h$','interpreter','latex')
title('Optimal Effort Levels','interpreter','latex')
legend('Benchmark','Lower $\kappa_v$','Location','northwest','interpreter','latex')
set(gca,'XTick',[], 'YTick', [])

subplot(1,2,2)
plot(c,S_e(N_z,:).*vcfunded(N_z,:),c,S_e_c(N_z,:).*vcfunded_c(N_z,:))
xlabel('Founder''s cost $c$','interpreter','latex')
ylabel('$S_e$','interpreter','latex')
title('Entrepreneur''s Surplus','interpreter','latex')
legend('Benchmark','Lower $\kappa_v$','Location','northwest','interpreter','latex')
set(gca,'XTick',[], 'YTick', [])
