function [val, grad, hess] = probit(b, Nk, T, X_vector, y_vector)

sign_y_all = 2 * y_vector - 1;

cvx_begin
        variable c(Nk, 1);


        cst = kron(c, ones(T,1) );
        Q =  - sum( log_normcdf( ( cst + X_vector * b ) .* sign_y_all ) );
        minimize( Q );
    cvx_end
    val = cvx_optval;
    c_out = c;



YL = reshape( y_vector, [T, Nk] );
YR = reshape( X_vector(:,1), [T Nk]);
x2 = reshape( X_vector(:,2), [T Nk]);

[logl, Grad, Hess,  Grho, Gbeta]=G_H(b, c_out,YL,YR,x2);

Hess( isnan(Hess)) = 0;
Hess( isinf(Hess)) = 0;
hess = -Hess;

Grad( isinf(Grad)) = 0;
Grad(isnan(Grad)) = 0;
grad =  - Grad;

end




