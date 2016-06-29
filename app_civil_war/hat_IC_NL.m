function [loglik, h_IC, group, post_c, post_b  ] = hat_IC_NL(N, y, X, group, a, b, K)
    T = size(y, 1)/N;
    ds = dataset( kron( (1:N)', ones(T, 1) ), repmat((1:T)', N, 1), y, X );
    ds.Properties.VarNames = {'N'  'T'  'y'  'X'};

global p
      %% post estimation
        post_c = zeros( N, 1 );
        post_b = zeros( size(b) );

        n_index = 1:N;
        QQ = zeros(K,1);
        for k = 1:K
            this_group = (group == k);
            this_group_index = n_index(this_group);
            this_N = length(this_group_index);
            n_extract = ismember(ds.N, this_group_index);
            
            this_y = ds.y(n_extract);
            this_X = ds.X(n_extract, :);
            sign_this_y = 2*this_y - 1;
                cvx_begin
                    variable c(this_N, 1);
                    variable b(1, p);

                    cst = kron(c, ones(T,1) );
                    Q =  -1/(N*T) * sum( log_normcdf( ( cst + this_X * b' ) .* sign_this_y ) );

                    minimize(Q);
                cvx_end
            post_b(this_group, :) = repmat(b, [this_N, 1] ) ;
            post_c(this_group) = c;
            QQ(k) = Q;
        end
        
        
    cst = kron(post_c, ones(T,1) ); 
    B = kron( post_b, ones(T,1)); % update the b as the post-estimator
    
    sign_y = 2*ds.y - 1;
    loglik = sum(QQ); 
    
    h_IC    = 0;
end
