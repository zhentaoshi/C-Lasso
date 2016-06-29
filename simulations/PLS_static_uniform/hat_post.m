function [ post_a, post_b ] = hat_post( y, X, b, a, group, K)

    N = size(b, 1);
    T = size(y, 1)/N;

    ds = dataset( kron( (1:N)', ones(T, 1) ), repmat((1:T)', N, 1), y, X );
    ds.Properties.VarNames = {'N'  'T'  'y'  'X'};


    post_a = zeros( size(a) );
    post_b = zeros( size(b) );

    for k = 1:K
        this_group = (group(:,k) == 1);
        post_a(k, :)  = post_est(N, ds, this_group);
        post_b(this_group, :) = repmat(post_a(k, :), sum( this_group ), 1) ;
    end

end
