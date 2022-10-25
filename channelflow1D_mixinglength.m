%% Fully developed channel flow solver using Prandtl'x mixing length model
% 46115 Turbuence modeling - Hamid Sarlak, Assoc. Prof. at DTU Wind and Energy Systems
% Disclaimer: this code is provided to you 
% without a guarantee on accuracy nor performance. 
% It is to be used as a basis to complete Assignment 2
% 
clc; clear; clf;
%% Channel Data
h = 0.1; %full channel height [m]
nu =1.48e-5; %kinematic viscosity [m^2/s]
rho = 1.225; % density [kg/m^3]

Re_tau = 180; % Re_tau = u_tau*h/nu, input paramerer 

%% Guess Theoretical Values

u_tau = Re_tau * nu /(h/2); %Re_tau = u_tau*h/nu, friction velocity [m/s]
dy_0_theory = 1/u_tau*nu; % y_plus = dy_0*u_tau/nu == 1, first grid spacing [m]
tau_w = u_tau^2*rho; %wall shear stress [Pa]
dp_dx_dimension = -2*tau_w/h; %[Pa/m]
dp_dx = dp_dx_dimension/(u_tau^2 * rho)*h; %[-] dimensionless pressure gradient

%% Grid Generation

R = 1.01; %stretching factor

dy_0 = dy_0_theory/(h/2);% non-dimensionalisation
N_y_viscous = log((R - 1)/dy_0 + 1)/log(R) + 1; %Calculate grid points required in viscous sublayer
N_y_viscous = ceil(N_y_viscous); %round to full value

viscous_lim = h/2; %limit for non-uniform grid (for h/2 the grid is always non-uniform)
y = [0, dy_0*cumsum(R.^(0:(N_y_viscous-2)))]; %no uniform grid in viscous sublayer [-]
y = y.*viscous_lim./max(y); %scaling of non-uniform grid to limits [m]

y = [y, linspace(y(end)+(y(end)-y(end-1)),h/2, (h/2 - y(end))/(y(end)-y(end-1)))]; %uniform grid from viscous_lim to h/2 [m] 
 
N_y = length(y); %total number of grid points [-]

% figure
% plot(0,y, 'k*'); %Mesh visualization
%% Solution Vector initialization

u = ones(1,N_y); %[m/s]
up = NaN(1,N_y); %[m/s]

%Boundary Condition
u(1) = 0; %[m/s]
%% Solving: Finite Difference

rms_err = 1; %initial value for rms_err
counter = 0; %counts number of iterations
tol_rms = 1e-5; %convergence criteria

kappa = 0.41; %van Karman constant

lm = NaN(1,N_y); %store for mixing length
VD = NaN(1,N_y); %store for Van Driest Damping function
nu_T = NaN(1,N_y); %store for kinemtic Viscosity

%Auxilary variables
A = NaN(1,N_y);
B = NaN(1,N_y);
C = NaN(1,N_y);

% Store stress profiles
viscous_stress = NaN(1,N_y);
turbulent_stress = NaN(1,N_y);

%Compute dimensionless quantities
l_plus = nu/u_tau; %[m]
y_plus = y./l_plus; %[-]

y = y./h; %[-], non dimensional y
u = u./u_tau; %[-], u_plus (non dimensional)

while  counter < 1e6  && rms_err > tol_rms
    
    for ii = 1:N_y
        
        if ii==1
            up(ii)=0;
            lm(ii) = 0;
            
        elseif ii < N_y
            
            %Van Driest Damping (Mixing length)
            A = 26;
            VD(ii) = (1-exp(-y_plus(ii)/A));
            lm(ii) = kappa * y(ii)*VD(ii);
            
            dudy = (u(ii+1) - u(ii-1))/(y(ii+1)-y(ii-1)); %Velocity Gradient
            
            nu_T(ii) = lm(ii)^2 * abs(dudy); %Eddy viscosity
            
            B(ii) = 1/(2*Re_tau) + 2*nu_T(ii); %auxialary variable
            C(ii) = (lm(ii)^2 - lm(ii-1)^2)/(y(ii)-y(ii-1))*dudy^2; %auxialary variable
            up(ii) = 1/(y(ii+1)-y(ii-1))*(u(ii+1)*(y(ii)-y(ii-1))+ u(ii-1)*(y(ii+1)-y(ii))...
                - 1/B(ii)*(dp_dx -C(ii))*0.5*(y(ii+1)-y(ii-1))*(y(ii+1)-y(ii))*(y(ii)-y(ii-1))); %calculate new velocity
            
            viscous_stress(ii) = nu*rho*dudy*u_tau/h;
            turbulent_stress(ii) = rho*nu_T(ii)*dudy*u_tau^2;
            
            
        else
            %up(ii) = up(ii-1); %forward differnce
            up(ii) = (-2*up(ii-1)+0.5*up(ii-1))*(-2/3); %second order
            lm(ii) = lm(ii-1);
        end
        
        
    end
    
    counter = counter +1;
    rms_err = sqrt(1/N_y.*sum((u-up).^2));
    u = up;
    
    %magnify to real solution (speed it up)
    u = u.*u_tau; % [m/s], use old u_tau
    tau_w = nu*rho * (u(2)-u(1))/(y(2)-y(1))*1/h; %[Pa]
    u_tau_new = sqrt(tau_w/rho); %[m/s], update
    u = u./u_tau_new; %[-]
   
end

%% Plotting
