function epsilon=dissipation(rho,P,mu,mu_T,y,sigma_epsilon,C_epsilon_1,C_epsilon_2,k,epsilon_old)
    N_y=length(y);

    A1=((2*mu+(mu_T(3:end)+mu_T(2:end-1))/sigma_epsilon)./(y(3:end)-y(2:end-1))./(y(3:end)-y(1:end-2)));
    B1=-((2*mu+(mu_T(3:end)+mu_T(2:end-1))/sigma_epsilon)./(y(3:end)-y(2:end-1))+(2*mu+(mu_T(2:end-1)+mu_T(1:end-2))/sigma_epsilon)./(y(2:end-1)-y(1:end-2)))./(y(3:end)-y(1:end-2));
    C1=((2*mu+(mu_T(2:end-1)+mu_T(1:end-2))/sigma_epsilon)./(y(2:end-1)-y(1:end-2))./(y(3:end)-y(1:end-2)));
    
    aWe=zeros(N_y,1);aPe=zeros(N_y,1);aEe=zeros(N_y,1);be=zeros(N_y,1);
    Ae=sparse(N_y,N_y);
    
    aEe(2:end-1)=A1;
    aPe(2:end-1)=B1;%+C_epsilon_1.*P./k(2:end-1)-2*C_epsilon_2*epsilon_old(2:end-1)./k(2:end-1);
    aWe(2:end-1)=C1;

    Ae=Ae+sparse(1,1,1,N_y,N_y,1);be(1)=0;
    Ae=Ae+sparse(N_y,N_y,-1.5,N_y,N_y,1);Ae=Ae+sparse(N_y,N_y-1,2,N_y,N_y,1);Ae=Ae+sparse(N_y,N_y-2,-0.5,N_y,N_y,1);be(N_y)=0;

    for i=2:N_y-1
        Ae=Ae+sparse(i,i-1,aWe(i),N_y,N_y,1);
        Ae=Ae+sparse(i,i,aPe(i),N_y,N_y,1);
        Ae=Ae+sparse(i,i+1,aEe(i),N_y,N_y,1);
        be(i)=rho*epsilon_old(i)/k(i)*(C_epsilon_2*epsilon_old(i)-C_epsilon_1*P(i-1));%C_epsilon_2/k(i)*epsilon_old(i)^2-C_epsilon_1/k(i)*P(i-1)*epsilon_old(i);
    end
    epsilon=Ae\be;
    
    
%     while error>tol
%         aEe(2:end-1)=A1;
%         aPe(2:end-1)=B1+C_epsilon_1.*P./k(2:end-1)-2*C_epsilon_2*epsilon_old(2:end-1)./k(2:end-1);
%         aWe(2:end-1)=C1;
%         
%         Ae=Ae+sparse(1,1,1,N_y,N_y,1);be(1)=epsilon1;
%         Ae=Ae+sparse(N_y,N_y,-1.5,N_y,N_y,1);Ae=Ae+sparse(N_y,N_y-1,2,N_y,N_y,1);Ae=Ae+sparse(N_y,N_y-2,-0.5,N_y,N_y,1);be(N_y)=0;
%         
%         for i=2:N_y-1
%             Ae=Ae+sparse(i,i-1,aWe(i),N_y,N_y,1);
%             Ae=Ae+sparse(i,i,aPe(i),N_y,N_y,1);
%             Ae=Ae+sparse(i,i+1,aEe(i),N_y,N_y,1);
%             be(i)=C_epsilon_2/k(i)*epsilon_old(i)^2-C_epsilon_1/k(i)*P(i-1)*epsilon_old(i);
%         end
%         epsilon=Ae\be;
%         error=norm(epsilon-epsilon_old)
%         epsilon_old=epsilon;
%     end
    
%     J=sparse(N_y-2,N_y-2);G=zeros(N_y-2,1);
%     error=1;tol = 10^-5;
%     
%     while error>tol
%         aEe(2:end-1)=A1;
%         aPe(2:end-1)=B1+C_epsilon_1.*P./k(2:end-1)-2*C_epsilon_2*epsilon(2:end-1)./k(2:end-1);
%         aWe(2:end-1)=C1;
%         
%         J=J+sparse(1,1,aPe(2),N_y-2,N_y-2,1);
%         J=J+sparse(1,2,aEe(2),N_y-2,N_y-2,1);
%         J=J+sparse(N_y-2,N_y-2,aPe(N_y-1),N_y-2,N_y-2,1);
%         J=J+sparse(N_y-2,N_y-3,aWe(N_y-1),N_y-2,N_y-2,1);
%         epsilon(end)=0.1;%(-2*epsilon(end-1)+0.5*epsilon(end-2))*(-2/3);
%         
%         G(1)=A1(1)*epsilon(3)+B1(1)*epsilon(2)+C1(1)*epsilon(1)+C_epsilon_1/k(2)*P(1)*epsilon(2)-C_epsilon_2/k(2)*epsilon(2)^2;
%         G(N_y-2)=A1(end)*epsilon(end)+B1(end)*epsilon(end-1)+C1(end)*epsilon(end-2)+C_epsilon_1/k(end-1)*P(end)*epsilon(end-1)-C_epsilon_2/k(end-1)*epsilon(end-1)^2;
%         
%         for i=3:N_y-2
%             J=J+sparse(i-1,i-2,aWe(i),N_y-2,N_y-2,1);
%             J=J+sparse(i-1,i-1,aPe(i),N_y-2,N_y-2,1);
%             J=J+sparse(i-1,i,aEe(i),N_y-2,N_y-2,1);
%             G(i-1)=A1(i-1)*epsilon(i+1)+B1(i-1)*epsilon(i)+C1(i-1)*epsilon(i-1)+C_epsilon_1/k(i)*P(i-1)*epsilon(i)-C_epsilon_2/k(i)*epsilon(i)^2;
%         end
%         delta=J\(-G);
%         epsilon(2:end-1)=epsilon(2:end-1)+delta;
%         error = norm(delta)
%     end
    
    
end