function omega=dissipation(rho,mu,mu_T,y,sigma_omega,alpha,beta,u_old,k_old,omega_old,F1)
    N_y=length(y);sigma_omega2=0.856;beta1=0.075;
    
    a1=(y(2:end-1)-y(1:end-2))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2));
    b1=((y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2))./(y(2:end-1)-y(1:end-2))-(y(2:end-1)-y(1:end-2))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2)));
    c1=-(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2))./(y(2:end-1)-y(1:end-2));
    
    dudy=u_old(3:end).*a1+u_old(2:end-1).*b1+u_old(1:end-2).*c1;
    dkdy=k_old(3:end).*a1+k_old(2:end-1).*b1+k_old(1:end-2).*c1;
    domegady=omega_old(3:end).*a1+omega_old(2:end-1).*b1+omega_old(1:end-2).*c1;

    A1=((2*mu+(mu_T(3:end).*sigma_omega(3:end)+mu_T(2:end-1).*sigma_omega(2:end-1)))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2)));
    B1=-((2*mu+(mu_T(3:end).*sigma_omega(3:end)+mu_T(2:end-1).*sigma_omega(2:end-1)))./(y(3:end)-y(2:end-1))+(2*mu+(mu_T(2:end-1).*sigma_omega(2:end-1)+mu_T(1:end-2).*sigma_omega(1:end-2)))./(y(2:end-1)-y(1:end-2)))./(y(3:end)-y(1:end-2));
    C1=((2*mu+(mu_T(2:end-1).*sigma_omega(2:end-1)+mu_T(1:end-2).*sigma_omega(1:end-2)))./(y(2:end-1)-y(1:end-2))./(y(3:end)-y(1:end-2)));
    
    aWe=zeros(N_y,1);aPe=zeros(N_y,1);aEe=zeros(N_y,1);be=zeros(N_y,1);
    Ae=sparse(N_y,N_y);
    
    aEe(2:end-1)=A1;
    aPe(2:end-1)=B1;
    aWe(2:end-1)=C1;

    Ae=Ae+sparse(1,1,1,N_y,N_y,1);be(1)=10*6*mu/rho/(beta1*(y(2)-y(1))^2);
    Ae=Ae+sparse(N_y,N_y,-1.5,N_y,N_y,1);Ae=Ae+sparse(N_y,N_y-1,2,N_y,N_y,1);Ae=Ae+sparse(N_y,N_y-2,-0.5,N_y,N_y,1);be(N_y)=0;

    for i=2:N_y-1
        Ae=Ae+sparse(i,i-1,aWe(i),N_y,N_y,1);
        Ae=Ae+sparse(i,i-0,aPe(i),N_y,N_y,1);
        Ae=Ae+sparse(i,i+1,aEe(i),N_y,N_y,1);
        be(i)=beta(i)*rho*omega_old(i)^2-alpha(i)*rho*dudy(i-1)^2-2*(1-F1(i))*sigma_omega2*rho/omega_old(i)*dkdy(i-1)*domegady(i-1);
    end
    omega=Ae\be;
    omega=max(omega,eps);

end