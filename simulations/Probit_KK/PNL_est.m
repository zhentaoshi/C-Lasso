
function [b_out, a_out, c_out] = PNL_est(N, T, b_initial, a_initial,  y, X, K, lam, R, tol )











global p


pen  = ones(N, K);
b_out  = repmat( b_initial, [ 1 1 K]) ;
a_out  = zeros(K, p);

b_old = ones(N, p);
a_old = zeros(1, p);

sign_y = (2*y - 1);

cvx_quiet(true)


if K == 1
    cvx_begin
    cvx_solver mosek

        variable c(N, 1);
        variable b(1, p);

        cst = kron(c, ones(T,1) );
        Q =  -1/(N*T) * sum( log_normcdf( ( cst + X * b' ) .* sign_y ) );

        minimize(Q);
    cvx_end

    a_out = b;
    b_out = b;
    c_out = c;

    else
    for r = 1:R
        for k = 1:K


            for kk = setdiff(1:K, k)
                pen(:, kk) =  norms( bsxfun(@minus, b_out(:, :, kk), a_out(kk, :) ), 2, 2 );
            end
            pen(:,k) = ones(N, 1);
            penalty_out  =  prod(pen, 2);

            cvx_begin
            cvx_solver mosek

                variable c(N, 1);
                variable b(N, p)
                variable a(1, p);

                B1 = kron( b, ones(T,1) );
                cst = kron(c, ones(T,1) );
                Q =  -1/(N*T) * sum( log_normcdf( ( cst + sum(X .* B1, 2) ) .* sign_y ) );


                pen_k = norms(  b - repmat(a, N, 1), 2, 2 );
                penalty =  penalty_out' * pen_k ;
                minimize( Q  + lam/N *  penalty);

            cvx_end

            pen(:,k) = pen_k;

            b_out(:, :, k) = b;
            a_out(k, :) = a;
            c_out = c;
        end

        a_new = a;
        b_new = b;

        if criterion( a_old, a_new, b_old, b_new, tol  ) == 1
            break;
        end


        a_old = a_new;
        b_old = b_new;

    end
end
end
