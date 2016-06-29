
global p N_cut a0



cvx_quiet true

N_cut =  N * [ 0.3, 0.6, 1 ];



group0 = zeros(N,1);
group0(1:N_cut(1)) = 1;
group0((N_cut(1)+1):N_cut(2)) = 2;
group0((N_cut(2)+1):N_cut(3)) = 3;


seed_r = seed;

proportion = [0.3, 0.3, 0.4];
a0_mat = repmat( a0(:,1)', [Rep, 1]);
POST_a1 = zeros(Rep, K); o_POST_a1 = POST_a1;
SD_a1 = zeros( Rep, K); o_SD_a1 = SD_a1;
SD_a1_rob = SD_a1; o_SD_a1_rob = SD_a1;
POST_a1_corr = POST_a1; o_POST_a1_corr = POST_a1;
HAT_a_corr_all = zeros(size(HAT_a_all));

tic
for r = 1:Rep
    r
    N = N_full;

    hat.group = HAT_group(:, r);

    seed_r = seed_r + 1;
    rng(seed_r);

    y = zeros(T, N);
    X = zeros(T, N, p);

    for kk = 1:3
        if kk == 1
            n_id = 1:N_cut(kk);
        else
            n_id = (N_cut(kk-1) + 1):N_cut(kk);
        end

        N_kk = length(n_id);
        [y(:, n_id), X(:, n_id, :) ] = DGP_NL_dyn(N_kk, T, a0(kk, :), a_constant(kk) );
    end



    unwanted = (mean(y) == 0) | (mean(y) ==1);
    wanted = ~unwanted;

    N = sum(wanted);





    alpha_hat0 = zeros(N,1);

    beta_hat0 = [ repmat( a0(1,:), [N_cut(1), 1] );...
        repmat( a0(2,:), [N_cut(2) - N_cut(1), 1] ); ...
        repmat( a0(3,:), [N_cut(3)-N_cut(2), 1] ) ];
    beta_hat0 = beta_hat0(wanted, :);




    for kk = 1:K




        this_group = logical (hat.group == kk);
        this.y = y(:, this_group);
        this.X = X(:, this_group, :);
        Nk = sum(this_group);

        y_vector = reshape(this.y, T*Nk, 1);
        X_vector = reshape(this.X, T*Nk, p);

        a_initial = HAT_a_all(kk, :, r)';


        [post_a, post_c] = solve(a_initial, Nk, T, X_vector, y_vector);

        SD_a1(r,kk) = SD_rho(post_a, post_c, Nk, T, this.y, this.X);
        SD_a1_rob(r,kk) = SD_rho_robust(post_a, post_c, Nk, T, this.y, this.X);
        HAT_a_corr_all(kk, :, r) = post_a;
        POST_a1(r, kk) = post_a(1);

        bias1 =  SPJ_zt( post_a, y_vector, X_vector, Nk, T)- post_a(1);

        POST_a1_corr(r,kk) = post_a(1) - bias1(1);
    end


    for kk = 1:K
        this_group = logical (group0 == kk);
        this.y = y(:, this_group);
        this.X = X(:, this_group, :);
        Nk = sum(this_group);

        y_vector = reshape(this.y, T*Nk, 1);
        X_vector = reshape(this.X, T*Nk, p);
        sign_y_all = 2*y_vector - 1;

        cvx_begin
        variable c(Nk, 1);
        variable b(p,1)

        B1 = kron( b, ones(T,1) );
        cst = kron(c, ones(T,1) );
        Q =  -1/(Nk*T) * sum( log_normcdf( ( cst + X_vector * b ) .* sign_y_all ) );
        minimize( Q );
        cvx_end

		o_initial  = a0(kk,:);
        [o_post_a, o_post_c] = solve(o_initial', Nk, T, X_vector, y_vector);

        o_SD_a1(r,kk) = SD_rho(o_post_a, o_post_c, Nk, T, this.y, this.X);
        o_SD_a1_rob(r,kk) = SD_rho_robust(o_post_a, o_post_c, Nk, T, this.y, this.X);
        o_POST_a1(r, kk) = o_post_a(1);

        o_bias1 =  SPJ_zt( o_post_a, y_vector, X_vector, Nk, T)-o_post_a(1);

        o_POST_a1_corr(r,kk) =  o_post_a(1) - o_bias1(1);

    end
    [POST_a1_corr(r,:), o_POST_a1_corr(r,:)];
    toc



    DEV_a1_corr(1:r,:) = POST_a1_corr(1:r,:) - a0_mat(1:r,:);
    o_DEV_a1_corr(1:r,:) =  o_POST_a1_corr(1:r,:) - a0_mat(1:r,:);



    mean( abs( DEV_a1_corr(1:r,:)./SD_a1(1:r,:)) <= 1.96 );
    mean( abs( o_DEV_a1_corr(1:r,:)./o_SD_a1(1:r,:)) <= 1.96 );


    bias = abs(mean(DEV_a1_corr(1:r,:) )) * proportion';
    RMSE = sqrt( mean( DEV_a1_corr(1:r,:).^2 ) * proportion');
    cover = mean( abs( DEV_a1_corr(1:r,:)./SD_a1_rob(1:r,:)) <= 1.96 ) * proportion';

    o_bias = abs( mean(o_DEV_a1_corr(1:r,:) ) ) * proportion';
    o_RMSE = sqrt( mean( o_DEV_a1_corr(1:r,:).^2 ) * proportion');
    o_cover = mean( abs( o_DEV_a1_corr(1:r,:)./o_SD_a1_rob(1:r,:)) <= 1.96 ) * proportion';

    summ = [RMSE, bias, cover, o_RMSE, o_bias, o_cover]


    if mod(r,min(50, Rep) ) == 0
        title = ['PNL_real_dyn_',num2str(N),'_T_', num2str(T),'_rep_', num2str(r),'_seed_', num2str(seed) ];
        export(dataset(summ), 'file', [title, '.csv'], 'Delimiter', ',');
        save( [title, '.mat'])
    end
end



DEV_a1_corr = POST_a1_corr - a0_mat;
o_DEV_a1_corr =  o_POST_a1_corr - a0_mat;

o_DEV_a1 = o_POST_a1 -a0_mat;
mean(o_DEV_a1)
mean(o_DEV_a1_corr)

mean( abs( DEV_a1_corr./SD_a1) <= 1.96 )
mean( abs( o_DEV_a1_corr./o_SD_a1) <= 1.96 )


bias = abs(mean(DEV_a1_corr )) * proportion';
RMSE = sqrt( mean( DEV_a1_corr.^2 ) * proportion');
cover = mean( abs( DEV_a1_corr./SD_a1_rob) <= 1.96 ) * proportion';



o_bias = abs( mean(o_DEV_a1_corr ) ) * proportion';
o_RMSE = sqrt( mean( o_DEV_a1_corr.^2 ) * proportion');
o_cover = mean( abs( o_DEV_a1_corr./o_SD_a1_rob) <= 1.96 ) * proportion';

summ = [RMSE, bias, cover, o_RMSE, o_bias, o_cover]


title = ['PNL_real_dyn_',num2str(N),'_T_', num2str(T),'_rep_', num2str(r),'_seed_', num2str(seed) ];
export(dataset(summ), 'file', [title, '.csv'], 'Delimiter', ',');
save( [title, '.mat'])
