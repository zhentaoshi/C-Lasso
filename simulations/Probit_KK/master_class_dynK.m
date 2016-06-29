global p N_cut a0 K_max

tic

N_full = N;
tol = 0.001;
R = 20;


p = 2;
N_cut =  N * [ 0.3, 0.6, 1 ];



group0 = zeros(N,1);
group0(1:N_cut(1)) = 1;
group0((N_cut(1)+1):N_cut(2)) = 2;
group0((N_cut(2)+1):N_cut(3)) = 3;




lambda1 = [0.0125,0.025,0.05, 0.1, 0.2];
numlam = length(lambda1);



a0 = [1, -1;
    .5, 0;
    0, 1];
a_constant = [-.5, -.25, 0];



rho = 1/4 * log(log(T)) / T;
IC_data = 2 * ones( K_max, numlam, Rep);
IC_dat = 2 * ones(Rep, K_max);

PPEE = zeros(numlam, Rep);
seed_r = seed;

for r = 1:Rep
    seed_r = seed_r + 1;
    rng(seed_r);
    N = N_full;
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
    beta_hat0 = beta_hat0(wanted,:);



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


    cvx_begin

    variable c(N, 1);
    variable b(p,1)

    B1 = kron( b, ones(T,1) );
    cst = kron(c, ones(T,1) );
    Q =  -1/(N*T) * sum( log_normcdf( ( cst + X * b ) .* (2*y-1) ) );
    minimize( Q );
    cvx_end

    IC_data(1,:,r) = 2 * cvx_optval;
    IC_dat(r,1) = 2 * cvx_optval;

    for ll = 1:numlam
        lam = lambda1(ll) *var(y)* T^(-1/3);

        for K = 2:K_max

            [b_K, a_out, c_out] = PNL_est(N, T, beta_hat0, alpha_hat0, y,  X, K, lam, R, tol );

            hat.a = a_out;
            hat.c = c_out;
            [hat.a, hat.b, hat.group] = report_b_coerce( b_K, hat.a, K);

            loglik = hat_IC_NL(N, y, X, hat.group, hat.b, hat.c, K);
            if K == 3;
                [PE] = mean( hat.group == group0(wanted) );
                PPEE( r ) = PE;        disp(PE)
            end


            IC_data(K, ll, r) =  2 * loglik;
        end
    end
    IC_dat(r,:) =  min(IC_data(:,:,r),[], 2);
    IC_final =  IC_dat(1:r,:) + rho * p * repmat( (1:K_max)-1, [r 1]);

    [ ~, hat_K] =  min(IC_final');
    disp( mean( hat_K == 3 ) );

    toc
    if mod(r,min(50, Rep) ) == 0
        title = ['IC_PNL_dyn_',num2str(N),'_T_', num2str(T),'_rep_', num2str(r),...
            '_seed_', num2str(seed) ];
        save(  [title, '.mat']);
    end
    r
end

PPEE;
[mean(PPEE,2), median(PPEE,2)]
title = ['IC_PNL_dyn_',num2str(N),'_T_', num2str(T),'_rep_', num2str(r),'_seed_', num2str(seed) ];
save(  [title, '.mat']);
