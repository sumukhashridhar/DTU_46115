function mu_T=eddy_viscosity(y,u_old,k_old,omega_old,rho,nu)
    a1=(y(2:end-1)-y(1:end-2))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2));
    b1=((y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2))./(y(2:end-1)-y(1:end-2))-(y(2:end-1)-y(1:end-2))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2)));
    c1=-(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2))./(y(2:end-1)-y(1:end-2));
    
    dudy=u_old(3:end).*a1+u_old(2:end-1).*b1+u_old(1:end-2).*c1;
    a1=0.31;beta_star=0.09;dudy=[0;dudy;0];
    arg2=max(2*sqrt(k_old)./(beta_star*omega_old.*y),500*nu./(y.^2.*omega_old));
    omega_vorticity=sqrt(2*(0.5*dudy).^2);
    F2=tanh(arg2.^2);
    mu_T=rho*a1*k_old./max(a1*omega_old,omega_vorticity.*F2);
    mu_T(1)=0;
end