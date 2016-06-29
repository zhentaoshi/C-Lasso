function [t] = PLS_super_coverage(N, T, Rep, seed)
rng(seed);
global p N_cut a0 R tol c_tol K0
cvx_solver mosek
p = 2;





N_cut =  N * [ 0.3, 0.6, 1 ];
K0 = length(N_cut);



a0 = [0.4, 1.6;      1, 1;    1.6, 0.4 ];



c_tol = 5;
tol = 0.0001;
R = 80;


sequence_lambda = 0.5 /( T^(1/3) );


myRep(Rep).H = 0;


beta_hat0 = zeros(N, p);
b0 = beta_hat0;


tic
cover_a = zeros(Rep,2);
cover_b = zeros(Rep,2);
cover_c = zeros(Rep,2);
A = zeros(3, Rep, 2);
B = zeros(N, Rep, 2);
C = zeros(N, Rep, 2);
II = zeros(Rep, 1);
for r = 1:Rep

    r
    toc
    [y, X] = DGP_static(N, T);


    for i = 1:N
        yy = reshape(y, [T N]);
        XX = reshape(X, [T N p]);
        b0(i,:) = a0( sum(N_cut < i) + 1, :);
        beta_hat0(i,:) = regress( yy(:, i) , permute( XX(:,i,:), [1 3 2]) );
    end



    Q = ones(1, 5) * 9999;

    lambda = sequence_lambda * var(y);
    [H, H_post] = outcome( N, T, beta_hat0, y, X, K0, lambda);
    e = y - sum(X.* kron( H_post.post_b, ones(T,1)) , 2);
    Q(3) = mean( e.^2 );

    myRep(r).H = H;
    myRep(r).H_post = H_post;
    A(:,r,1) = H.a(:,1);
    A(:,r,2) = H_post.post_a(:,1);
    B(:,r,1) = H.b_co(:,1);
    B(:,r,2) = H_post.post_b(:,1);


    single_a = abs( H.a(:,1) - a0(:,1) )./sqrt( H_post.vari_a(:,1) ) < 1.96;
    cover_a(r,1) =  [0.3 0.3 0.4] * single_a;
    single_a = abs( H_post.post_a(:,1) - a0(:,1) )./sqrt( H_post.vari_a_post(:,1) ) < 1.96;
    cover_a(r,2) =  [0.3 0.3 0.4] * single_a;

    single_b = abs( H.b_co(:,1) - b0(:,1) )./sqrt( H_post.vari_b(:,1) ) < 1.96;
    cover_b(r,1) = mean(single_b);
    single_b = abs( H_post.post_b(:,1) - b0(:,1) )./sqrt( H_post.vari_b_post(:,1) ) < 1.96;
    cover_b(r,2) = mean(single_b);




    for kk = [2 4 5]
        [~, H_post_K] = outcome_K( N, T, beta_hat0, y, X, kk, lambda);
        e = y - sum( X.*kron( H_post_K.post_b, ones(T,1)), 2);
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
        single_c = abs( H_K.b_co(:,1) - b0(:,1) )./sqrt( H_post.vari_b(:,1) ) < 1.96;
        cover_c(r,1) = mean(single_c);
        single_c = abs( H_post_K.post_b(:,1) - b0(:,1) )./sqrt( H_post_K.vari_b_post(:,1) ) < 1.96;
        cover_c(r,2) = mean(single_c);
    end
        [cover_a(r,:); cover_b(r,:); cover_c(r,:)]


end

RMSE(1,1) = [0.3 0.3 0.4] * mean( ( A(:,:,1) - repmat(a0(:,1), [1 Rep]) ).^2, 2);
RMSE(2,1) = [0.3 0.3 0.4] * mean( ( A(:,:,2) - repmat(a0(:,1), [1 Rep]) ).^2, 2);
RMSE(1,2) = mean( mean( ( B(:,:,1) - repmat(b0(:,1), [1 Rep]) ).^2, 2) );
RMSE(2,2) = mean( mean( ( B(:,:,2) - repmat(b0(:,1), [1 Rep]) ).^2, 2) );
RMSE(1,3) = mean( mean( ( C(:,:,1) - repmat(b0(:,1), [1 Rep]) ).^2, 2) );
RMSE(2,3) = mean( mean( ( C(:,:,2) - repmat(b0(:,1), [1 Rep]) ).^2, 2) );
RMSE = sqrt( RMSE );

BIAS(1,1) = [0.3 0.3 0.4] * (mean( A(:,:,1),2) - a0(:,1));
BIAS(2,1) = [0.3 0.3 0.4] * (mean( A(:,:,2),2) - a0(:,1));
BIAS(1,2) = mean( mean( B(:,:,1),2) - b0(:,1) ) ;
BIAS(2,2) = mean( mean( B(:,:,2),2) - b0(:,1) ) ;
BIAS(1,3) = mean( mean( C(:,:,1),2) - b0(:,1) ) ;
BIAS(2,3) = mean( mean( C(:,:,2),2) - b0(:,1) ) ;

out = [RMSE(:,1), BIAS(:,1), mean( cover_a )', ...
    RMSE(:,2), BIAS(:,2), mean( cover_b )', ...
    RMSE(:,3), BIAS(:,3), mean(cover_c)']

title = [ 'PLS_static_cover_N_', num2str(N), '_T_', num2str(T), '_Rep_', num2str(Rep),'_seed_', num2str(seed)];
save( [title, '.mat']);
export(dataset(out), 'file', [title, '.csv'], 'Delimiter', ',');

mean(II == 3)
t = toc
end
