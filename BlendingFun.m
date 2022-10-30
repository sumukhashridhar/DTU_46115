function F1=BlendingFun(y,k_old,omega_old,rho,nu)
    beta_star=0.09;sigma_omega2=0.856;
    
    a1=(y(2:end-1)-y(1:end-2))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2));
    b1=((y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2))./(y(2:end-1)-y(1:end-2))-(y(2:end-1)-y(1:end-2))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2)));
    c1=-(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2))./(y(2:end-1)-y(1:end-2));
    
    dkdy=k_old(3:end).*a1+k_old(2:end-1).*b1+k_old(1:end-2).*c1;
    domegady=omega_old(3:end).*a1+omega_old(2:end-1).*b1+omega_old(1:end-2).*c1;
    
    CD_kw=max(2*rho*sigma_omega2./omega_old(2:end-1).*dkdy.*domegady,1e-20);
    arg1=min(max(sqrt(k_old(2:end-1))./(beta_star*omega_old(2:end-1).*y(2:end-1)),500*nu./(y(2:end-1).^2.*omega_old(2:end-1))),4*rho*sigma_omega2.*k_old(2:end-1)./(CD_kw.*y(2:end-1).^2));
    F1=[1;tanh(arg1.^4);0];
end