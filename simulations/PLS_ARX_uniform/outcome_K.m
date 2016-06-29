function [H, H_post] = outcome_K( N, T, b0, y, X, K, lambda)










global  R tol

[b_K, a_unordered] = PLS_est_single_K(N, T, b0, y, X, K, lambda, R, tol);


[ H.b_co, H.group_co] = report_b( T, b_K, a_unordered, K);



TT = size(y, 1)/N;
ds = dataset( kron( (1:N)', ones(TT, 1) ), repmat((1:TT)', N, 1), y, X );
ds.Properties.VarNames = {'N'  'T'  'y'  'X'};

[H_post] = hat_IC_PLS(ds, H.b_co, a_unordered, K, H.group_co );

end
