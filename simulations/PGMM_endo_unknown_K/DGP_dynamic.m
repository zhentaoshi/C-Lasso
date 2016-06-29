function [y, X, Dy, DX, Z] = DGP_dynamic(N, T )




















    global N_cut a0 d p K0

    sigma_x = 1;
    sigma_e = [1 1 1];

    TT = T + d;




    y = zeros( TT, N );
    Dy = zeros( TT , N);
    Dy_lag = zeros(TT, N);
    Z = zeros( TT, N, d);

    if p > 1

        X_U = random('Normal', 0, sigma_x, [ TT+1  N  p-1 ] );
    end


    for i = 1:N
        kk = 1:K0;
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





                if t >= d + 2;
                    Dy    (t-1, i )   = y(t,i) - y(t-1, i);
                    Dy_lag(t-1, i )   = y( t - 1, i) - y( t - 2, i);
                    Z     (t-1, i, :) = y( (t - 2 - d + 1):(t - 2), i );
                end
            end
        end
    end



    DX = zeros( TT, N, p);
    DX(:, :, 1) = Dy_lag;
    if p > 1
        Lag_X_U = X_U(1:TT, :, :);
        D_X_U   = X_U(2:(TT+1), :, :) - Lag_X_U;
        DX(:, :, 2:p) = D_X_U;



        Z(:,:, (d+1):(p+d-1) ) =  D_X_U;

    end








    Dy(1:d, :) = [];
    DX(1:d, :, :) = [];
    Z(1:d, :, :) = [];



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




