function [t] = PGMM_super_coverage(N, T, Rep, seed)
rng(seed);
global p N_cut a0 R tol c_tol K0 d


p = 3;



title1 = [ 'PGMM_ARX_cover_N_', num2str(N), '_T_', num2str(T), '_Rep_', num2str(Rep),'_seed_', num2str(seed)];
load( [title1,'.mat'] );

d = 2;


N_cut =  N * [ 0.3, 0.6, 1 ];
K0 = length(N_cut);



a0 = [0.4, 1.6, 1.6;...
    0.6, 1, 1;...
    0.8, 0.4, 0.4 ];


a0_1 = a0(:,1);

a0_vector = zeros(N,1);
group0 = zeros(N, K0);
for k = 1:K0
    if k == 1
        group0(1: N_cut(k) , k) = 1;
        a0_vector(1:N_cut(k) ) = a0_1(k);
    else
        group0( (N_cut(k-1) + 1):N_cut(k), k ) = 1;
        a0_vector((N_cut(k-1) + 1):N_cut(k) ) = a0_1(k);
    end
end
group0 = logical(group0);



c_tol = 5;
tol = 0.0001;
R = 80;
cvx_solver mosek

sequence_lambda = .5 /( T^(1/3) );
beta_hat0 = zeros(N, p);
b0 = zeros(N,p);


tic
B = zeros(N, Rep, 4);
RMSE = zeros(Rep,4);
BIAS = zeros(Rep,4);

for r = 1:Rep

    r
    toc


    [~, ~, y, X, Z] = DGP_dynamic(N, T);


    for i = 1:N
        b0(i,:) = a0( sum(N_cut < i) + 1, :);
        Xt = permute( X(:,i,:), [1 3 2]);
        Zt = permute( Z(:,i,:), [1 3 2]);
        yt = y(:, i);
        beta_hat0(i,:) = pinv(Xt' * (Zt * Zt') * Xt) * (Xt' * (Zt * Zt') * yt);
    end



    lambda = sequence_lambda * var(reshape(y, 1, []) ) ;
    [H, H_post] = outcome_GMM( N, T, beta_hat0, y, X, Z, K0, lambda);
     b_PO = H_post.post_b(:,1);
     b_CL = H_post.b_corr(:,1);

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
    if I == 3
        B(:,r,3) = B(:,r,1);
        B(:,r,4) = B(:,r,2);
    else
        [H_K, H_post_K] = outcome_GMM_K( N, T, beta_hat0, y, X, Z, I, lambda);
        b_PO_K = H_post_K.post_b(:,1);
        b_CL_K = H_post_K.b_corr(:,1);
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
            a0_N_PO( b_PO_K== a_K_PO(ii)) = a0_perm_PO(ii);
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

title2 = [ 'PGMM_ARX_oracle2_N_', num2str(N), '_T_', num2str(T), '_Rep_', num2str(Rep),'_seed_', num2str(seed)];
save( [title2, '.mat']);
export(dataset(summ), 'file', [title, '.csv'], 'Delimiter', ',');

t = toc
end
