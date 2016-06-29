function [H] = hat_IC_PLS( ds, b, a, K, group)




      N = size(b, 1);
      T = size(ds,1)/N;
        post_a = zeros( size(a) );
        a_corr = post_a;
        post_a_corr = post_a;

        post_b = zeros( size(b) );
        b_corr = post_b;
        post_b_corr = post_b;

        for k = 1:K
            this_group = (group(:, k)  == 1);
            [a_corr(k, :), post_a(k, :), post_a_corr(k,:)]  = post_est_PLS(N, T, a(k,:), ds, this_group, 1);
            b_corr(this_group, :) = repmat(a_corr(k, :), sum( this_group ), 1) ;
            post_b(this_group, :) = repmat(post_a(k, :), sum( this_group ), 1) ;
            post_b_corr(this_group, :) = repmat(post_a_corr(k, :), sum( this_group ), 1) ;
        end
    H.a_corr = a_corr;
    H.post_a = post_a;
    H.post_a_corr = post_a_corr;
    H.b_corr = b_corr;
    H.post_b = post_b;
    H.post_b_corr = post_b_corr;
end
