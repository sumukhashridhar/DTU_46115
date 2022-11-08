clc
%clear
close all

%% Channel Data
h = 0.1; %full channel height [m]
nu =1.48e-5; %kinematic viscosity [m^2/s]
rho = 1.225; % density [kg/m^3]
mu=nu*rho;

Re_tau = 182; % Re_tau = u_tau*h/nu, input paramerer 
beta=0.05;%Relaxation factor, choose beta to be 0.05, 0.01 and 0.005 for Re=180 Re=590 Re=2000 respectively to avoid divergence
damp=0;%0 no damp. fun. 1 with damp. fun.

%% Guess Theoretical Values

u_tau = Re_tau * nu /h*2; %Re_tau = u_tau*h/nu, friction velocity [m/s]
dy_0_theory = 1/u_tau*nu; % y_plus = dy_0*u_tau/nu == 1, first grid spacing [m]
tau_w = u_tau^2*rho; %wall shear stress [Pa]
dp_dx = -tau_w*2/h; %[Pa/m]
% dp_dx = dp_dx_dimension/(u_tau^2 * rho)*h; %[-] dimensionless pressure gradient

%% Grid Generation

R = 1.01; %stretching factor

dy_0 = dy_0_theory/(h/2);% non-dimensionalisation
N_y_viscous = log((R - 1)/dy_0 + 1)/log(R) + 1; %Calculate grid points required in viscous sublayer
N_y_viscous = ceil(N_y_viscous); %round to full value

viscous_lim = h/2; %limit for non-uniform grid (for h/2 the grid is always non-uniform)
y = [0, dy_0*cumsum(R.^(0:(N_y_viscous-2)))]; %no uniform grid in viscous sublayer [-]
y = y.*viscous_lim./max(y); %scaling of non-uniform grid to limits [m]

y = [y, linspace(y(end)+(y(end)-y(end-1)),h/2, (h/2 - y(end))/(y(end)-y(end-1)))]'; %uniform grid from viscous_lim to h/2 [m] 
 
N_y = length(y); %total number of grid points [-]

% figure
% plot(0,y, 'k*'); %Mesh visualization
%% Solution Vector initialization
%Closure Coefficients
C_epsilon_1=1.44;C_epsilon_2=1.92;C_mu=0.09;sigma_k=1.0;sigma_epsilon=1.3;

% mu_T = zeros(N_y,1);
% u_old = mean_velocity(mu,mu_T,y,dp_dx); %[m/s]
[u_old,mu_T_old]=mixinglength(Re_tau,h,rho,nu);
epsilon_old=0.03*ones(N_y,1);
k_old=0.1*ones(N_y,1);

%Boundary Condition
% u(1) = 0; %[m/s]
%% Solving: Finite Difference

rms_err = 1; %initial value for rms_err
counter = 0; %counts number of iterations
tol_rms = 1e-10; %convergence criteria


% Store stress profiles
viscous_stress = NaN(1,N_y);
turbulent_stress = NaN(1,N_y);

%Compute dimensionless quantities
% l_plus = nu/u_tau; %[m]
% y_plus = y./l_plus; %[-]

% y = y./h; %[-], non dimensional y
% u = u./u_tau; %[-], u_plus (non dimensional)

aWr=zeros(N_y,1);aPr=zeros(N_y,1);aEr=zeros(N_y,1);br=zeros(N_y,1);
aWk=zeros(N_y,1);aPk=zeros(N_y,1);aEk=zeros(N_y,1);bk=zeros(N_y,1);
aWe=zeros(N_y,1);aPe=zeros(N_y,1);aEe=zeros(N_y,1);be=zeros(N_y,1);
Ar=sparse(N_y,N_y);Ak=sparse(N_y,N_y);Ae=sparse(N_y,N_y);

a1=(y(2:end-1)-y(1:end-2))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2));
b1=((y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2))./(y(2:end-1)-y(1:end-2))-(y(2:end-1)-y(1:end-2))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2)));
c1=-(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2))./(y(2:end-1)-y(1:end-2));

a2=2./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2));
b2=-2./(y(3:end)-y(1:end-2)).*(1./(y(3:end)-y(2:end-1))+1./(y(2:end-1)-y(1:end-2)));
c2=2./(y(2:end-1)-y(1:end-2))./(y(3:end)-y(1:end-2));

A1=((2*mu+mu_T_old(3:end)+mu_T_old(2:end-1))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2)));
B1=-((2*mu+mu_T_old(3:end)+mu_T_old(2:end-1))./(y(3:end)-y(2:end-1))+(2*mu+mu_T_old(2:end-1)+mu_T_old(1:end-2))./(y(2:end-1)-y(1:end-2)))./(y(3:end)-y(1:end-2));
C1=((2*mu+mu_T_old(2:end-1)+mu_T_old(1:end-2))./(y(2:end-1)-y(1:end-2))./(y(3:end)-y(1:end-2)));

acenter=((y(N_y)-y(N_y-1))/(y(N_y)-y(N_y-2))-(y(N_y)-y(N_y-2)/(y(N_y)-y(N_y-1))))/(y(N_y-2)-y(N_y-1));
bcenter=(y(N_y)-y(N_y-2))/(y(N_y)-y(N_y-1))/(y(N_y-2)-y(N_y-1));
ccenter=-(y(N_y)-y(N_y-1))/(y(N_y)-y(N_y-2))/(y(N_y-2)-y(N_y-1));

abc=2/(y(2)-y(3))*(1/(y(3)-y(1))-1/(y(2)-y(1)));
bbc=2/(y(2)-y(1))/(y(2)-y(3));
cbc=-2/(y(3)-y(1))/(y(2)-y(3));


y_plus=y*u_tau/nu;

while  counter < 1e6  && rms_err > tol_rms
    counter=counter+1;
    
    Rt=k_old.^2./nu./epsilon_old;Ry=sqrt(k_old).*y;
    f2=((1-0.3*exp(-Rt.^2))*damp+(1-damp)*1)*0+1;
    f_mu1=(1-exp(-y_plus/150)).^2*damp+(1-damp)*1;%exp(-3.4./(1+Rt/50).^2)*(damp)+(1-damp);
    
    mu_T=f_mu1.*eddy_viscosity(k_old,epsilon_old,rho,C_mu);
    
    u=mean_velocity(mu,mu_T,y,dp_dx);
    
    u=beta*u+(1-beta)*u_old;
    
    rms_err=norm(u-u_old)
    
    dudy=u(3:end).*a1+u(2:end-1).*b1+u(1:end-2).*c1;
    
    P=TK_production(mu_T,dudy,rho);
    
    k=TKE(rho,mu_T,mu,sigma_k,y,epsilon_old,P);
    
    k=beta*k+(1-beta)*k_old;
    
    epsilon=dissipation(rho,P,mu,mu_T,y,sigma_epsilon,C_epsilon_1,f2*C_epsilon_2,k,epsilon_old);
    
    epsilon=beta*epsilon+(1-beta)*epsilon_old;
    
    u_old=u;
    
    k_old=k;
    
    epsilon_old=epsilon;
    
end

%% Plotting
figure
plot(u,y,'LineWidth',2)
title('Mean velocity')
grid on

Re_stress=mu_T.*[0;dudy;0];
figure
plot(Re_stress,y,'LineWidth',2)
title('Reynolds stress')
xlabel('$\rho\overline{u''_iu''_j}$','Interpreter','latex')
ylabel('y')
grid on

figure
plot(mu_T,y)
title('Eddy viscosity \mu_T')
grid on

figure
plot(k,y,'LineWidth',2)
title('TKE')
xlim([0 Inf])
xlabel('$\frac{1}{2}\overline{u''_iu''_i}$','Interpreter','latex')
ylabel('y')
grid on

figure
plot(epsilon,y,'LineWidth',2)
title('TK dissipation rate \epsilon')
xlim([0 Inf])
xlabel('\epsilon')
ylabel('y')
grid on

u_plus=u/u_tau;
y1=linspace(y_plus(1),20,100);
figure
semilogx(y_plus,u_plus,'LineWidth',2)
hold on
semilogx(y1,y1,'--')
semilogx(y_plus,1/0.4187*log(9.793*y_plus),'--')
xlabel('y^+')
ylabel('U^+')
hold off
grid on