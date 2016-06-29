function [y, X] = DGP_static(N, T)












    global  p N_cut a0
    sigma_x = 1;
    sigma_e = 1;

    y = zeros(T, N );
    X = zeros(T, N, p);

    K = length(N_cut);

    for i = 1:N
        kk = 1:K;
        k = min( kk(i <= N_cut ) );

        XX = random('Normal', 0, sigma_x, T, p);
        e = random('Normal', 0, sigma_e, T, 1);
        yy = XX * a0(k, 1:p)' + e;

        X(:, i, :)  = permute( demean(XX), [1 3 2]) ;
        y(:, i) = demean( yy );
    end
    y = reshape(y, T*N, 1);
    X = reshape(X, T*N, p);
end

