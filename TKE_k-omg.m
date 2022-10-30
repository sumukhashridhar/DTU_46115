function k=TKE(rho,mu_T,mu,y,k,omega,P,sigma_k)
    N_y=length(y);
    beta_star=0.09;
    
    A1=((2*mu+(mu_T(3:end).*sigma_k(3:end)+mu_T(2:end-1).*sigma_k(2:end-1)))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2)));
    B1=-((2*mu+(mu_T(3:end).*sigma_k(3:end)+mu_T(2:end-1).*sigma_k(2:end-1)))./(y(3:end)-y(2:end-1))+(2*mu+(mu_T(2:end-1).*sigma_k(2:end-1)+mu_T(1:end-2).*sigma_k(1:end-2)))./(y(2:end-1)-y(1:end-2)))./(y(3:end)-y(1:end-2));
    C1=((2*mu+(mu_T(2:end-1).*sigma_k(2:end-1)+mu_T(1:end-2).*sigma_k(1:end-2)))./(y(2:end-1)-y(1:end-2))./(y(3:end)-y(1:end-2)));

    
    aWk=zeros(N_y,1);aPk=zeros(N_y,1);aEk=zeros(N_y,1);bk=zeros(N_y,1);
    Ak=sparse(N_y,N_y);

    aEk(2:end-1)=A1;
    aPk(2:end-1)=B1;
    aWk(2:end-1)=C1;
    Ak=Ak+sparse(1,1,1,N_y,N_y,1);bk(1)=0;
    Ak=Ak+sparse(N_y,N_y,-1.5,N_y,N_y,1);Ak=Ak+sparse(N_y,N_y-1,2,N_y,N_y,1);Ak=Ak+sparse(N_y,N_y-2,-0.5,N_y,N_y,1);bk(N_y)=0;
    for i=2:N_y-1
        Ak=Ak+sparse(i,i-1,aWk(i),N_y,N_y,1);
        Ak=Ak+sparse(i,i,aPk(i),N_y,N_y,1);
        Ak=Ak+sparse(i,i+1,aEk(i),N_y,N_y,1);
    end
    bk(2:end-1)=rho*(beta_star*omega(2:end-1).*k(2:end-1)-P);
    k=Ak\bk;
    k=max(k,eps);
end