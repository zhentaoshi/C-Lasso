function [beta_est, group, unclass] = report_b( b, a, K, tol )










    N = size(b, 1);

    group = zeros(N, K);
    est_group = zeros(N, K);
    for k = 1:K
        b_a = bsxfun(@minus, b(:, :, k), a(k, :) );
        group( ( sum( abs( b_a ), 2)  < tol ), k) = 1;
        est_group(:, k)  = group(:,k);
    end

    index_n = 1:N;
    in_group = index_n( sum(est_group, 2) == 1 );
    out_of_group = setdiff( (1:N), in_group);
    group = logical(group);
    unclass = N - sum( sum(est_group, 2) );




    if K == 1
        beta_est = b;
    else

        beta_est = 999 * ones( N, 2);

        for i = in_group
            for k = 1:K
                if group(i,k) == 1
                    beta_est(i,:) = a(k, :);
                end
            end
        end

        beta_est( out_of_group, : ) = mean( b( out_of_group, :, :), 3 );
    end
end
