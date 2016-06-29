function [Dy, DX, DZ] = DGP_endo(N, T )




















    global N_cut a0 d p K0 D
    D = p + d - 1;
    corr_coef = 0.3;



    TT = T + 1;


    y = zeros( TT, N );
    X1= zeros( TT, N );

    IV= random('Normal', 0, 1, [ TT N d] );
    Z2 = random('Normal', 0, 1, [TT  N ] );


    for i = 1:N
        kk = 1:K0;
        k = min( kk(i <= N_cut ) );

        mu = random('Normal', 0, 1, 1);
        IVi = permute(IV( :, i, : ), [ 1 3 2] );
        Z2i = Z2(:,  i) ;

        E = random('Normal', 0, 1, TT, 2);
        E = E * [1 corr_coef; corr_coef, 1];

        u = E(:, 1);
        e = E(:, 2);

        X1(:, i) = bsxfun(@plus, 0.2 * mu,   IVi * ones(2,1)/2 + 0.5*e) ;

        y(:,i) = mu + [ X1(:,i) Z2i] * a0(k , :)' + u;
    end

   X(:, :,1)  = X1;
   X(:, :,2)  = Z2;
   Z(:, :, 1:2) = IV;
   Z(:, :, 3) = Z2;



    Dy = y(2:TT, :) - y(1:T, :);
    DX = X(2:TT, :, :) - X(1:T, :, :);
    DZ = Z(2:TT, :, :) - Z(1:T, :, :);

end




