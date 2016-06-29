function [t] = PLS_super(N, T, Rep, seed)

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
    toc
     [y, X] = DGP_dynamic(N, T);
    TT = T + d;

    y = reshape(y, TT*N, 1);
    X = reshape(X, TT*N, p);


    dat = dataset( kron( (1:N)', ones(TT, 1) ), repmat((1:TT)', N, 1), y, X );
    dat.Properties.VarNames = {'N'  'T'  'y'  'X'};

    for k = 1:K
        this_group = (group0(:, k)  == 1);

        index = 1:N;
        g_index = index(this_group);
        g_data = dat( ismember(dat.N, g_index), : );
		thisY = g_data.y;
		thisX = g_data.X;

        post_a = regress( thisY, thisX );
        post_a_corr = post_a - bias_PLS(TT, post_a, thisY, thisX);

        post(k, :, r) = post_a_corr';
        post_var(k,:,r) = sqrt( diag( var_PLS(TT, post_a_corr, y, X)));
        cover(k,:,r) = abs( ( post(k, :, r) - a0(k,:) )./post_var(k,:,r) ) < 1.96;
    end
end

myTitle = [ 'PLS_ARX_Oracle_', num2str(N), '_T_', num2str(T), '_Rep_', num2str(Rep),'_seed_', num2str(seed), '.mat'];
save( myTitle);

t = toc;
end
