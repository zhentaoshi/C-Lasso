function [H, H_post] = outcome_GMM( N, T, b0, y, X, Z, K, lambda)










global a0 R tol K0 p d

[b_K, a_unordered] = PGMM_est(N, T, b0, y, X, Z, K0, lambda, R, tol);



perm_hat = perm_a( a0, a_unordered);
H.a = a_unordered(perm_hat, :);
b_K = b_K(:,:, perm_hat);

[~, H.b_co, ~, H.group_co] = report_b( T, b_K, H.a, K);


TT = size(y,1);

        X = reshape(X,[TT*N p]);
        Z = reshape(Z,[TT*N p + d - 1]);
        y = reshape(y,[TT*N 1]);



ds = dataset( kron( (1:N)', ones(TT, 1) ), repmat((1:TT)', N, 1), y, X, Z );
ds.Properties.VarNames = {'N'  'T'  'y'  'X' 'Z'};

[H_post] = hat_IC_PGMM(ds, H.b_co, H.a, K, H.group_co );

end
