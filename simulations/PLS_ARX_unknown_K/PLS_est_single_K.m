function [b_out, a_out] = PLS_est_single_K(N, T, b0, y, X, K, lambda, R, tol)


    global p

    pen  = ones(N, K);
    b_out  = repmat(b0, [ 1 1 K]);
    a_out  = zeros(K, p);

    b_old = ones(N, p);
    a_old = zeros(1, p);

    TT = T + 3;
    cvx_quiet(true)

        for r = 1:R
            for k = 1:K


                for kk = setdiff(1:K, k)
                    pen(:, kk) =  norms( bsxfun(@minus, b_out(:, :, kk), a_out(kk, :) ), 2, 2 );
                end
                pen(:,k) = ones(N, 1);
                penalty_out  =  prod(pen, 2);

                cvx_begin
                    variable b(N, p)
                    variable a(1, p);
                    B1 = kron( b, ones(TT,1) );
                    Q =  1/(N*TT) * sum_square(  y - sum(X .* B1, 2) );


                    pen_k = norms(  b - repmat(a, N, 1), 2, 2 );
                    penalty =  penalty_out' * pen_k ;


                    minimize( Q  + lambda/N *  penalty);
                cvx_end
                pen(:,k) = pen_k;

                b_out(:, :, k) = b;
                a_out(k, :) = a;

            end



            b_new = b_out(:, :, k);
            a_new = a_out(k, :);


            if criterion( a_old, a_new, b_old, b_new, tol  ) == 1
                break;
            end


            a_old = a_out(k, :);
            b_old = b_out(:, :, k);

        end
    end
