function mu_T = eddy_viscosity(lm,dudy,rho)
    mu_T = rho*lm.^2.*abs(dudy);
    mu_T = [0;mu_T;(-2*mu_T(end-1)+0.5*mu_T(end-2))*(-2/3)];
end