K0 = 3;
a0_1 = a0(:,1);


B1_true =   [ repmat( a0(1,:), [N_cut(1), 1] );...
    repmat( a0(2,:), [N_cut(2) - N_cut(1), 1] ); ...
    repmat( a0(3,:), [N_cut(3)-N_cut(2), 1] ) ];
B1_known_CL  = repmat( B1_true(:,1), [1 Rep]);
B1_known_PO = B1_known_CL;

BIAS1 = POST_a1 - POST_a1_corr;
BIAS1 = BIAS1';


B = zeros(N, Rep, 4);
RMSE = zeros(Rep,4);
BIAS = zeros(Rep,4);


for r = 1:Rep
    for kk = 1:3
        a_CL = HAT_a(kk, r)-BIAS1(kk,r);
        B1_known_CL( HAT_group(:,r ) == kk, r) = a_CL;
        a_PO = POST_a1_corr(r, kk);
        B1_known_PO( HAT_group(:,r ) == kk, r) = a_PO;
    end
end

B1_unknown_CL = B1_known_CL;
B1_unknown_PO = B1_known_PO;

lambda1 = 0.05;


seed_r = seed;
for r = 1:Rep

    seed_r = seed_r + 1;
    rng(seed_r);
    N = N_full;


   b_CL = B1_known_CL( :, r );
   b_PO = B1_known_PO( :, r );


        a_CL = unique( b_CL );
        group = zeros(N, K0);
        for ii = 1:K0
            group( b_CL== a_CL(ii), ii ) = 1;

        end

        [~, order] = sort( sum(group), 2, 'descend' );


        a_CL =  a_CL( order );

        a0_perm_CL = perm_a_K(a_CL, a0_1);

        a_PO = unique( b_PO );
        a_PO =  a_PO( order );
        a0_perm_PO = perm_a_K(a_PO, a0_1);

        a0_N_CL = zeros(N, 1);
        a0_N_PO = a0_N_CL;

        for ii = 1:K0
            a0_N_CL( b_CL == a_CL(ii) ) = a0_perm_CL(ii);
            a0_N_PO( b_PO== a_PO(ii)) = a0_perm_PO(ii);
        end


        B(:,r,1) = b_PO - a0_N_PO;
        B(:,r,2) = b_CL - a0_N_CL;


    I = II(r);
    if I== 3
        B(:,r,3) = B(:,r,1);
        B(:,r,4) = B(:,r,2);
    else
        r

        K = II(r);
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
        y_vector = reshape(y, T*N, 1);
        X_vector = reshape(X, T*N, p);

        lam = lambda1 *var(y_vector)* T^(-1/3);


        [b_K, a_out, c_out] = PNL_est(N, T, beta_hat0, alpha_hat0, y_vector,  X_vector, K, lam, R, tol );
        a_out(a_out >= 10 ) = 10;
        a_out(a_out <=-10 ) = -10;
        hat.a = a_out;
        hat.c = c_out;
        [hat.a, hat.b, hat.group] = report_b_coerce( b_K, hat.a, K);




    for kk = 1:K




        this_group = logical (hat.group == kk);
        this.y = y(:, this_group);
        this.X = X(:, this_group, :);
        Nk = sum(this_group);

        y_vector = reshape(this.y, T*Nk, 1);
        X_vector = reshape(this.X, T*Nk, p);

        a_initial = hat.a(kk, :)';


        [post_a, post_c] = solve(a_initial, Nk, T, X_vector, y_vector);
        bias1 =  SPJ_zt( post_a, y_vector, X_vector, Nk, T)- post_a(1);

        B1_unknown_CL( hat.group == kk , r) = hat.a(kk,1) - bias1(1);
        B1_unknown_PO( hat.group == kk , r) = post_a(1) - bias1(1);
    end


    b_PO_K = B1_unknown_PO( : , r);
    b_CL_K = B1_unknown_CL( : , r);


        a_K_CL = unique( b_CL_K );
        group_K = zeros(N, I);
        for ii = 1:I
            group_K( b_CL_K== a_K_CL(ii), ii ) = 1;

        end

        [~, order] = sort( sum(group_K), 2, 'descend' );


        a_K_CL =  a_K_CL( order );

        a0_perm_CL = perm_a_K(a_K_CL, a0_1);

        a_K_PO = unique( b_PO_K );
        a_K_PO =  a_K_PO( order );
        a0_perm_PO = perm_a_K(a_K_PO, a0_1);

        a0_N_CL = zeros(N, 1);
        a0_N_PO = a0_N_CL;

        for ii = 1:I
            a0_N_CL( b_CL_K == a_K_CL(ii) ) = a0_perm_CL(ii);
            a0_N_PO( b_PO_K == a_K_PO(ii)) = a0_perm_PO(ii);
        end

        B(:,r,3) = b_PO_K - a0_N_PO;
        B(:,r,4) = b_CL_K - a0_N_CL;
    end

	for jj = 1:4
		RMSE(r,jj) = sqrt( mean( mean(  B(:,1:r,jj).^2) ) ) ;
		BIAS(r,jj) = mean( mean( B(:,1:r,jj))  ) ;
		summ=[  RMSE(r,:), BIAS(r,:)];
	end
	summ

end

bs = permute( abs( mean(B) ), [2 3 1]);
summ =  [RMSE(Rep, :), mean(bs) ]


title = ['PNL_dyn_unknown_',num2str(N),'_T_', num2str(T),'_rep_', num2str(r),'_seed_', num2str(seed) ];
save(  [title, '.mat']);
export( dataset(summ), 'file', [title,'.csv'], 'Delimiter', ',')
