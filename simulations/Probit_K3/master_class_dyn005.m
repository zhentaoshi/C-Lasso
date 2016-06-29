global p N_cut a0

tic

p = 2;
N_full = N;
HAT_group = zeros(N, Rep);
HAT_a = zeros(3,Rep);
HAT_a_all = zeros(3,p, Rep);








N_cut =  N * [ 0.3, 0.6, 1 ];



group0 = zeros(N,1);
group0(1:N_cut(1)) = 1;
group0((N_cut(1)+1):N_cut(2)) = 2;
group0((N_cut(2)+1):N_cut(3)) = 3;

numlam = 1;
lamb.grid = numlam;
lamb.min  = 0.02;
lamb.max  = 0.20;
lambda1 = lamb.min * (lamb.max / lamb.min ).^( ( (1:lamb.grid) - 1) /( lamb.grid -1 ) );
lambda1 = 0.05;


a0 = [1, -.5;
     .5, 0;
    0.2, .5];
a_constant = [0, .25, .5];




tol = 0.001;
R = 20;

PPEE = zeros(numlam, Rep);

seed_r = seed;
for r = 1:Rep
    N = N_full;

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
    y = y(:, wanted);
    X = X(:, wanted, :);

    alpha_hat0 = zeros(N,1);

    beta_hat0 = [ repmat( a0(1,:), [N_cut(1), 1] );...
        repmat( a0(2,:), [N_cut(2) - N_cut(1), 1] ); ...
        repmat( a0(3,:), [N_cut(3)-N_cut(2), 1] ) ];
    beta_hat0 = beta_hat0(wanted, :);



    X_all= reshape(X, T*N, p);

    y_all = reshape(y, T*N,1);
    beta_temp = fitglm( X_all , y_all, 'distr', 'binomial', 'link', 'probit' );
    beta_temp_coef = beta_temp.Coefficients.Estimate;
    y_pred = ([ X_all, ones(size(X_all,1),1)] * beta_temp_coef >= 0);
    pseudo_R = mean(y_pred .* y_all);


    for i = 1:N
        yi = y(:,i);
        if 0.2 <= mean(yi) && mean(yi) <= 0.8
            xi = permute( X(:,i,:), [1 3 2]);
            beta_temp = fitglm( xi, yi, 'distr', 'binomial', 'link', 'probit' );
            beta_temp_coef = beta_temp.Coefficients.Estimate;

            if max(abs(beta_temp_coef)) <= 5
                alpha_hat0(i) = beta_temp_coef(1);
                beta_temp_coef(1) = [];
                beta_hat0(i,:) = beta_temp_coef;
            end
        end
    end
    y = reshape(y, T*N, 1);
    X = reshape(X, T*N, p);

    lambda = lambda1 *var(y)* T^(-1/3);
    for lam = lambda

        for K = 3

            [b_K, a_out, c_out] = PNL_est(N, T, beta_hat0, alpha_hat0, y,  X, K, lam, R, tol );
            a_out(a_out >= 10 ) = 10;
            a_out(a_out <=-10 ) = -10;
            hat.a = a_out;
            hat.c = c_out;
            [hat.a, hat.b, hat.group] = report_b_coerce( b_K, hat.a, K);


            PE = mean( hat.group == group0(wanted) )
			HAT_group( wanted,  r) = hat.group;
            HAT_a(:,r) = hat.a(:,1);
            HAT_a_all(:, :, r) = hat.a;
            PPEE(lam == lambda, r ) = PE

        end

    end

    toc

    r
end

PPEE
[mean(PPEE,2), median(PPEE,2)]
N = N_full;
title = ['PNL_K3_dyn_',num2str(N),'_T_', num2str(T),'_rep_', num2str(r),'_seed_', num2str(seed) ];
save(  [title, '.mat'] );
