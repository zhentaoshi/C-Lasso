function [hat] = PLS( N, T, y, X, lambda, rho, K, R, tol)


[b_K, hat.a] = PLS_est_single_K(N, T, y, X, K, lambda, R, tol);
[hat.b, hat.group, U] = report_b( b_K, hat.a, K, tol);
[~, hat.IC, hat.post_a, hat.post_b] = hat_IC(y, X, hat.b, hat.a, U, K, rho);


[hat.Ek, hat.Fk, hat.Nk] = report(hat.group, hat.post_a, K);
end

