clc
clear
close all

%% Channel Data
h = 0.1; %full channel height [m]
nu =1.48e-5; %kinematic viscosity [m^2/s]
rho = 1.225; % density [kg/m^3]
mu=nu*rho;

Re_tau = 182; % Re_tau = u_tau*h/nu, input paramerer 

%% Guess Theoretical Values

u_tau = Re_tau * nu /h; %Re_tau = u_tau*h/nu, friction velocity [m/s]
dy_0_theory = 1/u_tau*nu; % y_plus = dy_0*u_tau/nu == 1, first grid spacing [m]
tau_w = u_tau^2*rho; %wall shear stress [Pa]
dp_dx = -tau_w/h; %[Pa/m]
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

mu_T = zeros(N_y,1);
u_old = mean_velocity(mu,mu_T,y,dp_dx); %[m/s]
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

A1=((2*mu+mu_T(3:end)+mu_T(2:end-1))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2)));
B1=-((2*mu+mu_T(3:end)+mu_T(2:end-1))./(y(3:end)-y(2:end-1))+(2*mu+mu_T(2:end-1)+mu_T(1:end-2))./(y(2:end-1)-y(1:end-2)))./(y(3:end)-y(1:end-2));
C1=((2*mu+mu_T(2:end-1)+mu_T(1:end-2))./(y(2:end-1)-y(1:end-2))./(y(3:end)-y(1:end-2)));

acenter=((y(N_y)-y(N_y-1))/(y(N_y)-y(N_y-2))-(y(N_y)-y(N_y-2)/(y(N_y)-y(N_y-1))))/(y(N_y-2)-y(N_y-1));
bcenter=(y(N_y)-y(N_y-2))/(y(N_y)-y(N_y-1))/(y(N_y-2)-y(N_y-1));
ccenter=-(y(N_y)-y(N_y-1))/(y(N_y)-y(N_y-2))/(y(N_y-2)-y(N_y-1));

abc=2/(y(2)-y(3))*(1/(y(3)-y(1))-1/(y(2)-y(1)));
bbc=2/(y(2)-y(1))/(y(2)-y(3));
cbc=-2/(y(3)-y(1))/(y(2)-y(3));

alpha=0.1;beta=0.05;

while  counter < 1e6  && rms_err > tol_rms
    counter=counter+1;
    
    mu_T=eddy_viscosity(k_old,epsilon_old,rho,C_mu);
    
    u=mean_velocity(mu,mu_T,y,dp_dx);
    
    u=beta*u+(1-beta)*u_old;
    
    rms_err=norm(u-u_old)
    
    dudy=u(3:end).*a1+u(2:end-1).*b1+u(1:end-2).*c1;
    
    P=TK_production(mu_T,dudy,rho);
    
    k=TKE(rho,mu_T,mu,sigma_k,y,epsilon_old,P);
    
    k=beta*k+(1-beta)*k_old;
    
%     epsilon1=0;%2*mu/rho*(abc*k(1)+bbc*k(2)+cbc*k(3));
    
    epsilon=dissipation(rho,P,mu,mu_T,y,sigma_epsilon,C_epsilon_1,C_epsilon_2,k,epsilon_old);
    
    epsilon=beta*epsilon+(1-beta)*epsilon_old;
    
    u_old=u;
    
    k_old=k;
    
    epsilon_old=epsilon;
    
end

%% Plotting
figure
plot(u,y)
title('Mean velocity')

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
plot(epsilon,y)
title('TK dissipation rate \epsilon')