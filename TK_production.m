function P=TK_production(mu_T,dudy,rho)
    P=mu_T(2:end-1)/rho.*dudy.^2;
end