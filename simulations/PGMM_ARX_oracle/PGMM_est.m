function [b_out, a_out] = PGMM_est(N, T, b0, y, X, Z, K, lambda, R, tol)















    global p d



    pen  = ones(N, K);
    b_out  = repmat( b0, [1 1 K] );
    a_out  = zeros(K, p);

    b_old = ones(N, p);
    a_old = zeros(1, p);

    cvx_quiet(true)


        for r = 1:R
            for k = 1:K


                for kk = setdiff(1:K, k)
                    pen(:, kk) =  norms( bsxfun(@minus, b_out(:, :, kk), a_out(kk, :) ), 2, 2 );
                end
                pen(:,k) = ones(N, 1);
                penalty_out  =  prod(pen, 2);

                cvx_begin
                cvx_solver sedumi
                    variable b(N, p)
                    variable a(1, p);

                    if p == 1
                        B1 = kron( b', ones(T, 1) );
                        e = y - X .* B1;
                    else
                        B1 = repmat( permute( b, [3 1 2] ), [T 1 1] );
                        e = y - sum( X .* B1, 3);
                    end
                    g = repmat(e, [1 1 p+d-1 ] ) .* Z ;
                    g = mean( g, 1 );
                    g = reshape( g, N*(p-1+d) , 1 );
                    Q = (1/N) * sum_square(g);


                    pen_k = norms(  b - repmat(a, N, 1), 2, 2 );
                    penalty =  penalty_out' * pen_k ;


                    minimize( Q  + lambda/N *  penalty);
                cvx_end
                pen(:,k) = pen_k;

                b_out(:, :, k) = b;
                a_out(k, :) = a;
            end




            a_new = a;
            b_new = b;

            a_numerator = sum( sum( (a_new - a_old).^2) );
            a_denominator = sum( sum( a_old.^2) ) + 0.0001;

            b_numerator = sum(sum( sum( (b_new - b_old).^2)) );
            b_denominator = sum( sum( sum( b_old.^2) )) + 0.0001;

            dist1=  a_numerator / a_denominator;
            dist2=  b_numerator / b_denominator;

            if (dist1< tol && dist2< tol )
                break
            end


            a_old = a_new;
            b_old = b_new;

        end
end
