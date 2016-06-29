function [beta_est_co, group_co] = report_b( T, b, a, K )










global p tune_tol c_tol


    N = size(b, 1);
    distance = zeros(N,K);

    tune_tol = c_tol/( sqrt(N*T) * log( log(N*T) ));

    group_co = zeros(N, K);
    for k = 1:K

        nominator = norms( bsxfun(@minus, b(:, :, k), a(k, :) ), 2 ,2 );


        distance(:,k) = nominator./norm( a(k,:) ) ;

    end
    [~, which] = min(distance, [], 2);





    beta_est_co = 999 * ones(N, p );

    for i = 1:N
        group_co(i, which(i)) = 1;
        beta_est_co(i, :)  = a(which(i), :);
    end


end
