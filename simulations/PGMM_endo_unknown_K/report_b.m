function [beta_est, beta_est_co, group, group_co, unclass] = report_b( T, b, a, K )










global p tune_tol c_tol


    N = size(b, 1);
    distance = zeros(N,K);

    tune_tol = c_tol/( sqrt(N*T) * log( log(N*T) ));

    group = zeros(N, K);
    group_co = zeros(N, K);
    for k = 1:K

        nominator = norms( bsxfun(@minus, b(:, :, k), a(k, :) ), 2 ,2 );


        distance(:,k) = nominator./norm( a(k,:) ) ;

    end
    [min_val, which] = min(distance, [], 2);




    beta_est = 999 * ones(N, p );
    beta_est_co = 999 * ones(N, p );

    for i = 1:N
        group_co(i, which(i)) = 1;
        beta_est_co(i, :)  = a(which(i), :);
        if  min_val(i) < tune_tol
            group(i, which(i) ) = 1;
        end
    end


    index_n = 1:N;
    in_group = index_n( sum(group, 2) == 1 );
    out_of_group = setdiff( (1:N), in_group );
    unclass = length(out_of_group);




    for i = in_group
        for k = 1:K
            if group(i,k) == 1
                beta_est(i,:) = a(k, :);
            end
        end
    end

    beta_est( out_of_group, : ) = mean( b( out_of_group, :, :), 3 );
end
