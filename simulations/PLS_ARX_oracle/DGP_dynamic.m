function [y, X] = DGP_dynamic(N, T )




















    global N_cut a0 d p K

    sigma_x = 1;
    sigma_e = [1 1 1];

    TT = T + d;




    y = zeros( TT, N );
    Z = zeros( TT, N, d);


        X_U = random('Normal', 0, sigma_x, [ TT+1  N  p-1 ] );


    for i = 1:N
        kk = 1:K;
        k = min( kk(i <= N_cut ) );
        e = random('Normal', 0, 1, TT+1, 1);
        for t = 1: (TT+1)
            mu = random('Normal', 0, 1, 1);
            c = mu*(1-a0(k,1));

            X_U_it = bsxfun(@plus, 0.2*mu, X_U(t, i, :) ) ;

            X_U_it = permute(X_U_it, [3 1 2]);



            if t == 1
                y(t, i ) = mu + a0(k, 2:p) * X_U_it + sigma_e(k) * e(t);

            elseif t > 1
                y(t, i ) = c + a0(k, 1:p) * [y(t - 1, i); X_U_it ] + sigma_e(k) * e(t);


            end
        end
    end


    X = zeros( TT + 1, N, p);
    X(2:(TT+1), :, 1 ) = y( 1:TT,: );
    if p > 1
        X(2:(TT+1) , :, 2:p) = X_U(2:(TT+1), :, :) ;
    end

    y(  1, : ) = [];
    X(  1, :, :) = [];


    y = demean(y);
    for  pp = 1:p
        X(:, :, pp) = demean( X(:, :, pp) );
    end
end




