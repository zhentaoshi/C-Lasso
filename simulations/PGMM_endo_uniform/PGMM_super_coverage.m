function [t] = PGMM_super_coverage(N, T, Rep, seed)
rng(seed);
global p N_cut a0 R tol c_tol K0 d
p = 2;



d = 2;


N_cut =  N * [ 0.3, 0.6, 1 ];
K0 = length(N_cut);



a0 = [ .2, 1.8;...
     1, 1;...
     1.8, 0.2 ];



c_tol = 5;
tol = 0.0001;
R = 80;

cvx_solver mosek

sequence_lambda = .5 /( T^(1/3) );
beta_hat0 = zeros(N, p);
b0 = zeros(N,p);


tic
tic
cover_a = zeros(Rep,4);
cover_b = zeros(Rep,4);
cover_c = zeros(Rep,4);
A = zeros(3, Rep, 4);
B = zeros(N, Rep, 4);
C = zeros(N, Rep, 4);
II = zeros(Rep, 1);

for r = 1:Rep

    r
    toc
    Q = ones(1, 5) * 9999;

    [y, X, Z] = DGP_endo(N, T);


    for i = 1:N
        b0(i,:) = a0( sum(N_cut < i) + 1, :);
        Xt = permute( X(:,i,:), [1 3 2]);
        Zt = permute( Z(:,i,:), [1 3 2]);
        yt = y(:, i);
        beta_hat0(i,:) = pinv(Xt' * (Zt * Zt') * Xt) * (Xt' * (Zt * Zt') * yt);
    end



    lambda = sequence_lambda * var(reshape(y, 1, []) ) ;
    [H, H_post] = outcome_GMM( N, T, beta_hat0, y, X, Z, K0, lambda);

    yy = reshape(y, [N*T 1]);
    XX = reshape(X, [N*T,p]);
    e  = yy - sum(XX.* kron( H_post.post_b, ones(T,1)) , 2);

    Q(3) = mean( e.^2 ) ;


    A(:,r,1) = H.a(:,1);
    A(:,r,2) = H_post.post_a(:,1);
    A(:,r,3) = H_post.a_corr(:,1);
    A(:,r,4) = H_post.post_a_corr(:,1);

    B(:,r,1) = H.b_co(:,1);
    B(:,r,2) = H_post.post_b(:,1);
    B(:,r,3) = H_post.b_corr(:,1);
    B(:,r,4) = H_post.post_b_corr(:,1);



    single_a = abs( H.a(:,1) - a0(:,1) )./sqrt( H_post.vari_a(:,1) ) < 1.96;
    cover_a(r,1) =  [0.3 0.3 0.4] * single_a;
    single_a = abs( H_post.post_a(:,1) - a0(:,1) )./sqrt( H_post.vari_a_post(:,1) ) < 1.96;
    cover_a(r,2) =  [0.3 0.3 0.4] * single_a;
    single_a = abs( H_post.a_corr(:,1) - a0(:,1) )./sqrt( H_post.vari_a_post(:,1) ) < 1.96;
    cover_a(r,3) =  [0.3 0.3 0.4] * single_a;
    single_a = abs( H_post.post_a_corr(:,1) - a0(:,1) )./sqrt( H_post.vari_a_post(:,1) ) < 1.96;
    cover_a(r,4) =  [0.3 0.3 0.4] * single_a;

    single_b = abs( H.b_co(:,1) - b0(:,1) )./sqrt( H_post.vari_b(:,1) ) < 1.96;
    cover_b(r,1) = mean(single_b);
    single_b = abs( H_post.post_b(:,1) - b0(:,1) )./sqrt( H_post.vari_b_post(:,1) ) < 1.96;
    cover_b(r,2) = mean(single_b);
    single_b = abs( H_post.b_corr(:,1) - b0(:,1) )./sqrt( H_post.vari_b_post(:,1) ) < 1.96;
    cover_b(r,3) = mean(single_b);





    for kk = [2 4 5]
        [~, H_post_K] = outcome_GMM_K( N, T, beta_hat0, y, X, Z, kk, lambda);
        e = yy - sum( XX.*kron( H_post_K.post_b, ones(T,1)), 2);
        Q(kk) =  mean(e.^2)  ;


    end
    [~, I ] = min( log(Q) + 2/3 * (N*T)^(-.5) * p * (1:5) )
    II(r) = I;
    if I == 3
        C(:,r,:) = B(:,r,:);
        cover_c(r,:) = cover_b(r,:);
    else
        [H_K, H_post_K] = outcome_GMM_K( N, T, beta_hat0, y, X, Z, I, lambda);
        C(:,r,1) = H_K.b_co(:,1);
        C(:,r,2) = H_post_K.post_b(:,1);
        C(:,r,3) = H_post_K.b_corr(:,1);
        C(:,r,4) = H_post_K.post_b_corr(:,1);

        single_c = abs( H_K.b_co(:,1) - b0(:,1) )./sqrt( H_post.vari_b(:,1) ) < 1.96;
        cover_c(r,1) = mean(single_c);
        single_c = abs( H_post_K.post_b(:,1) - b0(:,1) )./sqrt( H_post_K.vari_b_post(:,1) ) < 1.96;
        cover_c(r,2) = mean(single_c);
        single_c = abs( H_post_K.b_corr(:,1) - b0(:,1) )./sqrt( H_post_K.vari_b_post(:,1) ) < 1.96;
        cover_c(r,3) = mean(single_c);


    end
    [cover_a(r,:); cover_b(r,:); cover_c(r,:)]
end


RMSE = zeros(3, 4);
BIAS = zeros(3, 4);

for jj = 1:4
    RMSE(jj,1) = [0.3 0.3 0.4] * mean( ( A(:,:,jj) - repmat(a0(:,1), [1 Rep]) ).^2, 2);
    RMSE(jj,2) = mean( mean( ( B(:,:,jj) - repmat(b0(:,1), [1 Rep]) ).^2, 2) );
    RMSE(jj,3) = mean( mean( ( C(:,:,jj) - repmat(b0(:,1), [1 Rep]) ).^2, 2) );
    BIAS(jj,1) = [0.3 0.3 0.4] * (mean( A(:,:,jj),2) - a0(:,1));
    BIAS(jj,2) = mean( mean( B(:,:,jj),2) - b0(:,1) ) ;
    BIAS(jj,3) = mean( mean( C(:,:,jj),2) - b0(:,1) ) ;
end

RMSE = sqrt( RMSE );


out = [RMSE(:,1), BIAS(:,1), mean( cover_a )', ...
    RMSE(:,2), BIAS(:,2), mean( cover_b )', ...
    RMSE(:,3), BIAS(:,3), mean(cover_c)'];

disp(out)

mean(II == 3)

title = [ 'PGMM_endo_cover_N_', num2str(N), '_T_', num2str(T), '_Rep_', num2str(Rep),'_seed_', num2str(seed)];
save( [title, '.mat']);
export(dataset(out), 'file', [title, '.csv'], 'Delimiter', ',');

t = toc
end
