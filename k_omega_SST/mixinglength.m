function [u,mu_T]=mixinglength(Re_tau,h,rho,nu)
    mu = nu*rho;
    
    u_tau = Re_tau * nu /(h/2); %Re_tau = u_tau*h/nu, friction velocity [m/s]
    dy_0_theory = 1/u_tau*nu; % y_plus = dy_0*u_tau/nu == 1, first grid spacing [m]
    tau_w = u_tau^2*rho; %wall shear stress [Pa]
    dp_dx = -2*tau_w/h; %[Pa/m]
    
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

mu_T_old=u_old;
%% Solving: Finite Difference

rms_err = 1; %initial value for rms_err
counter = 0; %counts number of iterations
tol_rms = 1e-5; %convergence criteria

kappa = 0.41; %van Karman constant


%Compute dimensionless quantities
l_plus = nu/u_tau; %[m]
y_plus = y./l_plus; %[-]

% y = y./h; %[-], non dimensional y
% u = u./u_tau; %[-], u_plus (non dimensional)

a1=(y(2:end-1)-y(1:end-2))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2));
b1=((y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2))./(y(2:end-1)-y(1:end-2))-(y(2:end-1)-y(1:end-2))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2)));
c1=-(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2))./(y(2:end-1)-y(1:end-2));

X_old=zeros(2*N_y,1);

X_old(1:2:end)=u_old;X_old(2:2:end)=mu_T_old;

while  counter < 1e6  && rms_err > tol_rms     
            
    A1=zeros(N_y,1);B1=zeros(N_y,1);C1=zeros(N_y,1);
    A1(2:end-1)=((2*mu+mu_T_old(3:end)+mu_T_old(2:end-1))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2)));
    B1(2:end-1)=-((2*mu+mu_T_old(3:end)+mu_T_old(2:end-1))./(y(3:end)-y(2:end-1))+(2*mu+mu_T_old(2:end-1)+mu_T_old(1:end-2))./(y(2:end-1)-y(1:end-2)))./(y(3:end)-y(1:end-2));
    C1(2:end-1)=((2*mu+mu_T_old(2:end-1)+mu_T_old(1:end-2))./(y(2:end-1)-y(1:end-2))./(y(3:end)-y(1:end-2)));
    
    
    J=sparse(2*N_y,2*N_y);G=zeros(2*N_y,1);
            
    %Van Driest Damping (Mixing length)
    A = 26;
    VD = (1-exp(-y_plus/A));
    lm = kappa * y.*VD;
    
    dudy=u_old(3:end).*a1+u_old(2:end-1).*b1+u_old(1:end-2).*c1; %Velocity Gradient
    
    J=J+sparse(1,1,1,2*N_y,2*N_y,1);G(1)=u_old(1)-0;
    J=J+sparse(2,2,1,2*N_y,2*N_y,1);G(2)=mu_T_old(1)-0;
    
    for ii=2:N_y-1
        J=J+sparse(2*(ii-1)+1,2*(ii-1)-1,C1(ii),2*N_y,2*N_y,1);
        J=J+sparse(2*(ii-1)+1,2*(ii-1)+1,B1(ii),2*N_y,2*N_y,1);
        J=J+sparse(2*(ii-1)+1,2*(ii-1)+3,A1(ii),2*N_y,2*N_y,1);
        G(2*(ii-1)+1)=A1(ii)*u_old(ii+1)+B1(ii)*u_old(ii)+C1(ii)*u_old(ii-1)-dp_dx;
        
        J=J+sparse(2*(ii-1)+2,2*(ii-1)+2,1     ,2*N_y,2*N_y,1);
        G(2*(ii-1)+2)=mu_T_old(ii)-rho*lm(ii).^2 .* abs(dudy(ii-1));
    end
    
    for ii=2*N_y-1:2*N_y
        J=J+sparse(ii,ii-4,-0.5,2*N_y,2*N_y,1);
        J=J+sparse(ii,ii-2, 2.0,2*N_y,2*N_y,1);
        J=J+sparse(ii,ii-0,-1.5,2*N_y,2*N_y,1);
    end
    G(2*N_y-1)=-0.5*u_old(N_y-2)+2*u_old(N_y-1)-1.5*u_old(N_y);
    G(2*N_y-0)=-0.5*mu_T_old(N_y-2)+2*mu_T_old(N_y-1)-1.5*mu_T_old(N_y);
    
    delta=J\(-G);
    
    X=X_old+0.1*delta;
    
    counter = counter +1;
    rms_err = norm(delta);
%     u_old = u;
    X_old=X;
    u_old=X(1:2:end);
    mu_T_old=X(2:2:end);
end

u=u_old;
mu_T=mu_T_old;

end