function u=mean_velocity(mu,mu_T,y,dp_dx)
    N_y=length(y);

    A1=((2*mu+mu_T(3:end)+mu_T(2:end-1))./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2)));
    B1=-((2*mu+mu_T(3:end)+mu_T(2:end-1))./(y(3:end)-y(2:end-1))+(2*mu+mu_T(2:end-1)+mu_T(1:end-2))./(y(2:end-1)-y(1:end-2)))./(y(3:end)-y(1:end-2));
    C1=((2*mu+mu_T(2:end-1)+mu_T(1:end-2))./(y(2:end-1)-y(1:end-2))./(y(3:end)-y(1:end-2)));

    aWr=zeros(N_y,1);aPr=zeros(N_y,1);aEr=zeros(N_y,1);br=zeros(N_y,1);
    Ar=sparse(N_y,N_y);

    aEr(2:end-1)=A1;
    aPr(2:end-1)=B1;
    aWr(2:end-1)=C1;
    
    % Compute mean velocity
    Ar=Ar+sparse(1,1,1,N_y,N_y,1);br(1)=0;
    Ar=Ar+sparse(N_y,N_y,-1.5,N_y,N_y,1);Ar=Ar+sparse(N_y,N_y-1,2,N_y,N_y,1);Ar=Ar+sparse(N_y,N_y-2,-0.5,N_y,N_y,1);br(N_y)=0;
    for i=2:N_y-1
        Ar=Ar+sparse(i,i-1,aWr(i),N_y,N_y,1);
        Ar=Ar+sparse(i,i,aPr(i),N_y,N_y,1);
        Ar=Ar+sparse(i,i+1,aEr(i),N_y,N_y,1);
        br(i)=dp_dx;
    end
    u=Ar\br;
end