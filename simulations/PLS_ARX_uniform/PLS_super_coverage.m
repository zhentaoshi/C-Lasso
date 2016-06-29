function [t] = PLS_super_coverage(N, T, Rep, seed)
rng(seed);
global p N_cut a0 R tol c_tol K0 d
p = 3;



d = 3;


N_cut =  N * [ 0.3, 0.6, 1 ];
K0 = length(N_cut);



a0 = [0.4, 1.6, 1.6;...
    0.6, 1, 1;...
    0.8, 0.4, 0.4 ];



c_tol = 5;
tol = 0.0001;
R = 80;
cvx_solver mosek

sequence_lambda = .5 /( T^(1/3) );


beta_hat0 = zeros(N, p);
b0 = beta_hat0;


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
    [y, X] = DGP_dynamic(N, T);


    for i = 1:N
        b0(i,:) = a0( sum(N_cut < i) + 1, :);
        beta_hat0(i,:) = regress( y(:, i) , permute( X(:,i,:), [1 3 2]) );
    end

    y = reshape(y, (T+d)*N, 1);

    X = reshape(X, (T+d)*N, p);


    Q = ones(1, 5) * 9999;

    lambda = sequence_lambda * var(y);
    [H, H_post] = outcome( N, T, beta_hat0, y, X, K0, lambda);
    e = y - sum(X.* kron( H_post.post_b, ones(T+3,1)) , 2);
    Q(3) = mean( e.^2 );

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
    single_b = abs( H_post.post_b_corr(:,1) - b0(:,1) )./sqrt( H_post.vari_b_post(:,1) ) < 1.96;
    cover_b(r,4) = mean(single_b);




    for kk = [2 4 5]
        [~, H_post_K] = outcome_K( N, T, beta_hat0, y, X, kk, lambda);
        e = y - sum( X.*kron( H_post_K.post_b, ones(T+3,1)), 2);
        Q(kk) = mean(e.^2) ;
    end
    [~, I ] = min( log(Q) + 2/3 * (N*T)^(-.5) * p * (1:5) )
    II(r) = I;
    if I == 3
        C(:,r,:) = B(:,r,:);
        cover_c(r,:) = cover_b(r,:);
    else
        [H_K, H_post_K] = outcome_K( N, T, beta_hat0, y, X, I, lambda);
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
        single_c = abs( H_post_K.post_b_corr(:,1) - b0(:,1) )./sqrt( H_post_K.vari_b_post(:,1) ) < 1.96;
        cover_c(r,4) = mean(single_c);
    end
        [cover_a(r,:); cover_b(r,:); cover_c(r,:)]


end

RMSE = zeros(4, 3);
BIAS = zeros(4, 3);

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

title = [ 'PLS_ARX_cover_N_', num2str(N), '_T_', num2str(T), '_Rep_', num2str(Rep),'_seed_', num2str(seed)];
save( [title, '.mat']);
export(dataset(out), 'file', [title, '.csv'], 'Delimiter', ',');

mean(II == 3)
t = toc
end
