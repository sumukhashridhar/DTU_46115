clc
clear
close all

%% Channel Data
h = 0.1; %full channel height [m]
nu =1.48e-5; %kinematic viscosity [m^2/s]
rho = 1.225; % density [kg/m^3]
mu=nu*rho;

Re_tau = 180; % Re_tau = u_tau*h/nu, input paramerer 
Beta=0.005;%Relaxation factor, choose beta to be 0.05, 0.01 and 0.005 for Re=180 Re=590 Re=2000 respectively to avoid divergence

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
sigma_k1=0.85;sigma_k2=1;
sigma_omega1=0.5;sigma_omega2=0.856;
alpha1=5/9;alpha2=0.44;
beta1=0.075;beta2=0.0828;

% mu_T = zeros(N_y,1);
% u_old = mean_velocity(mu,mu_T,y,dp_dx); %[m/s]
[u_old,mu_T_old]=mixinglength(Re_tau,h,rho,nu);
omega_old=0.3*ones(N_y,1);
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

damp=1;%0 no damp. fun. 1 with damp. fun.
y_plus=y*u_tau/nu;

%% Solving finite difference equation
while  counter<1e10 && rms_err>tol_rms
    counter=counter+1;
    
    F1=BlendingFun(y,k_old,omega_old,rho,nu);
    sigma_k=F1*sigma_k1+(1-F1)*sigma_k2;
    sigma_omega=F1*sigma_omega1+(1-F1)*sigma_omega2;
    alpha=F1*alpha1+(1-F1)*alpha2;
    beta=F1*beta1+(1-F1)*beta2;
    
    u=mean_velocity(mu,mu_T_old,y,dp_dx);
    
    u=Beta*u+(1-Beta)*u_old;
    
    rms_err=norm(u-u_old)

    P=TK_production(mu_T_old,y,u_old,rho,k_old,omega_old);
    
    k=TKE(rho,mu_T_old,mu,y,k_old,omega_old,P,sigma_k);
    
    k=Beta*k+(1-Beta)*k_old;
    
    omega=dissipation(rho,mu,mu_T_old,y,sigma_omega,alpha,beta,u_old,k_old,omega_old,F1);
    
    omega=Beta*omega+(1-Beta)*omega_old;
    
    mu_T=eddy_viscosity(y,u,k,omega,rho,nu);
    
    u_old=u;
    k_old=k;
    omega_old=omega;
    mu_T_old=mu_T;
    
    
    
end

%% Plotting
figure
plot(u,y)
title('Mean velocity')

dudy=u_old(3:end).*a1+u_old(2:end-1).*b1+u_old(1:end-2).*c1;
Re_stress=mu_T.*[0;dudy;0];
figure
plot(Re_stress,y)
title('Reynolds stress')

figure
plot(mu_T,y)
title('Eddy viscosity \mu_T')

figure
plot(k,y)
title('TKE')

figure
plot(omega,y)
title('TK dissipation rate \epsilon')

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