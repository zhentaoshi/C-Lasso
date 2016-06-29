function [v] = IC( sigma, rho, p,  K)
    v = log(sigma) + p * rho' * K;
end
