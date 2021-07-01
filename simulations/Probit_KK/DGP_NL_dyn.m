function [Y, X] = DGP_NL_dyn(N, T, a0, const)




    a0 = a0';
    TT = 100+T;

    sigma_x = 1;

    y = zeros(N,TT);
    ystar = y;

    e =  randn(N,TT);
    mu = randn(N,1);

    x2 = repmat( 0.1*mu, [1, TT]) + random('Normal', 0, sigma_x, N, TT);


    y0 = binornd(1, 0.5, N, 1);
    for t = 1:TT
        if t == 1
            ystar(:,t) = [y0, x2(:, t) ] * a0 + 0.1*mu+ const + e(:,t);
            y(:,t) = (ystar(:,t) > 0);
        else
            ystar(:,t) = [y(:, t-1), x2(:, t) ] * a0 + 0.1*mu+ const + e(:,t);
            y(:,t) = (ystar(:,t) > 0);
        end
    end

    Y = y( :, (TT - T + 1):TT );
    X(:,:,1) = y(:, (TT - T):(TT-1));
    X(:,:,2) = x2(:, (TT - T + 1):TT);

    Y = Y';
    X = permute(X, [2 1 3]);

end

