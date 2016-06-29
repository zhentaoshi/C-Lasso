function [v] = IC_NL( loglik, rho, p,  K)
    v = loglik + rho.*p*K;
end
