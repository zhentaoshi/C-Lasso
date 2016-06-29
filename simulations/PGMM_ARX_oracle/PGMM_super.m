function [t] = PGMM_super(N, T, Rep, seed)

rng(seed);
global p N_cut a0 K d
p = 3;
d = 3;
N_cut =  N * [ 0.3, 0.6,  1 ];
K = length(N_cut);





a0 = [0.4, 1.6, 1.6;...
    0.6, 1, 1;...
    0.8, 0.4, 0.4 ];

group0 = zeros(N, K);
for k = 1:K
    if k == 1
        group0(1: N_cut(k) , k) = 1;
    else
        group0( (N_cut(k-1) + 1):N_cut(k), k ) = 1;
    end
end
group0 = logical(group0);


post = zeros(K, p , Rep);
post_var = post;
cover   = post;

tic
for r = 1:Rep

    r
    toc;
    [~, ~, y, X, Z] = DGP_dynamic(N, T, K);


    X = reshape(X,[T*N p]);
    Z = reshape(Z,[T*N p + d - 1]);
    y = reshape(y,[T*N 1]);

    dat = dataset( kron( (1:N)', ones(T, 1) ), repmat((1:T)', N, 1), y, X, Z );
    dat.Properties.VarNames = {'N'  'T'  'y'  'X' 'Z'};

    for k = 1:K
        this_group = (group0(:, k)  == 1);

        index = 1:N;
        g_index = index(this_group);
        g_data = dat( ismember(dat.N, g_index), : );

        y = g_data.y;
        X = g_data.X;
        Z = g_data.Z;

        XZZX = (X' * Z) * (Z' * X);
        XZZy = (X' * Z) * (Z' * y);
        post_a = (pinv(XZZX) * XZZy)';

        post(k, :, r) = post_a;
        post_var(k,:,r) = sqrt( diag( var_post_GMM(T, post_a', y, X, Z)));
        cover(k,:,r) = abs( ( post(k, :, r) - a0(k,:) )./post_var(k,:,r) ) < 1.96;
    end
end

myTitle = [ 'PGMM_ARX_Oracle_', num2str(N), '_T_', num2str(T), '_Rep_', num2str(Rep),'_seed_', num2str(seed), '.mat'];
save( myTitle);

t = toc;
end
