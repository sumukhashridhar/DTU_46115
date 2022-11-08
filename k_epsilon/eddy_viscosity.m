function mu_T=eddy_viscosity(k,epsilon,rho,C_mu)
    mu_T=rho*C_mu*k.^2./epsilon;
end