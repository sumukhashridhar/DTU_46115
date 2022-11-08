%% Fully developed channel flow solver using Prandtl'x mixing length model
% 46115 Turbuence modeling - Hamid Sarlak, Assoc. Prof. at DTU Wind and Energy Systems
% Disclaimer: this code is provided to you 
% without a guarantee on accuracy nor performance. 
% It is to be used as a basis to complete Assignment 2
% 
close all
clc; 
%clear;
%% Channel Data
h = 0.1; %full channel height [m]
nu =1.48e-5; %kinematic viscosity [m^2/s]
rho = 1.225; % density [kg/m^3]
mu = nu*rho;

Re_tau = 8000; % Re_tau = u_tau*h/nu, input paramerer 

%% Guess Theoretical Values

u_tau = Re_tau * nu /(h/2); %Re_tau = u_tau*h/nu, friction velocity [m/s]
dy_0_theory = 1/u_tau*nu; % y_plus = dy_0*u_tau/nu == 1, first grid spacing [m]
tau_w = u_tau^2*rho; %wall shear stress [Pa]
dp_dx = -2*tau_w/h; %[Pa/m]
% dp_dx = dp_dx_dimension/(u_tau^2 * rho)*h; %[-] dimensionless pressure gradient

%% Grid Generation

R = 1.01; %stretching factor

dy_0 = dy_0_theory/(h/2);% non-dimensionalisation
N_y_viscous = log((R - 1)/dy_0 + 1)/log(R) + 1; %Calculate grid points required in viscous sublayer
N_y_viscous = ceil(N_y_viscous); %round to full value

viscous_lim = h/2; %limit for non-uniform grid (for h/2 the grid is always non-uniform)
y = [0, dy_0*cumsum(R.^(0:(N_y_viscous-2)))]; %no uniform grid in viscous sublayer [-]
y = y.*viscous_lim./max(y); %scaling of non-uniform grid to limits [m]

y = [y, linspace(y(end)+(y(end)-y(end-1)),h/2, (h/2 - y(end))/(y(end)-y(end-1)))]; %uniform grid from viscous_lim to h/2 [m] 

y=y';
N_y = length(y); %total number of grid points [-]

% figure
% plot(0,y, 'k*'); %Mesh visualization
%% Solution Vector initialization

u_old = ones(N_y,1); %[m/s]

%Boundary Condition
u_old(1) = 0; %[m/s]
%% Solving: Finite Difference

rms_err = 1; %initial value for rms_err
counter = 0; %counts number of iterations
tol_rms = 1e-8; %convergence criteria

kappa = 0.41; %van Karman constant

lm = zeros(N_y,1); %store for mixing length
VD = NaN(N_y,1); %store for Van Driest Damping function
nu_T = NaN(N_y,1); %store for kinemtic Viscosity


% Store stress profiles
viscous_stress = zeros(1,N_y);
turbulent_stress = zeros(1,N_y);

%Compute dimensionless quantities
l_plus = nu/u_tau; %[m]
y_plus = y./l_plus; %[-]

% y = y./h; %[-], non dimensional y
% u = u./u_tau; %[-], u_plus (non dimensional)

a1=(y(2:end-1)-y(1:end-2))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2));
b1=((y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2))./(y(2:end-1)-y(1:end-2))-(y(2:end-1)-y(1:end-2))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2)));
c1=-(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2))./(y(2:end-1)-y(1:end-2));

while  counter < 1e6  && rms_err > tol_rms     
            
    %Van Driest Damping (Mixing length)
    A = 26;
    VD = (1-exp(-y_plus/A));
    lm = kappa * y.*VD;
    
    dudy=u_old(3:end).*a1+u_old(2:end-1).*b1+u_old(1:end-2).*c1; %Velocity Gradient
    
    nu_T = lm(2:end-1).^2 .* abs(dudy); %Eddy viscosity
    nu_T=[0;nu_T;0];
    
    mu_T = rho*nu_T;
    
    u=mean_velocity(mu,mu_T,y,dp_dx);
    
%     viscous_stress = nu*rho*dudy;
%     turbulent_stress = rho*nu_T.*dudy;


    
    counter = counter +1;
    rms_err = rms(u_old-u)
    u_old = u;
    
   
end

%% Plotting
figure
plot(u,y)
title('Mean velocity')

figure
plot(mu_T,y,mu_T*0+mu,y)
legend('\mu_T','\mu')
title('Eddy viscosity \mu_T')

u_plus=u/u_tau;
y1=linspace(y_plus(1),20,100);

figure
semilogx(y_plus,u_plus,'LineWidth',2)
hold on
semilogx(y1,y1,'--', 'HandleVisibility','off')
semilogx(y_plus,1/0.4187*log(9.793*y_plus),'--', 'HandleVisibility','off')
xlabel('y^+')
ylabel('U^+')
hold off
grid on
