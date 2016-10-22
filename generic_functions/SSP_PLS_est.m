function [b_est, a_out, group_est] = SSP_PLS_est(N, T, y, X, K, lambda, R)
% Su, Shi and Phillips (2017)

% PLS estimation by the iterative algorithm

% INPUT:
%   N
%   T
%   y: (TN * 1)
%   X: (TN * p)
%   K: number of groups to be classified
%   R: maximum number of iterations
%   tol: tolerence level to judge convergence

% OUTPUT:
%   b_est: estimated beta  (N*p)
%   b_out: estiamted alpha (K*p)
%   group_est: estimated group identity

p = size(X,2);
tol = 0.0001; % tolerence level to judge convergence

% use the within-group regression as the initial value
beta_hat0 = zeros(N, p);
for i = 1:N
    yy = reshape(y, [T N]);
    XX = reshape(X, [T N p]);
    beta_hat0(i,:) = regress( yy(:, i) , permute( XX(:,i,:), [1 3 2]) );
end


pen  = ones(N, K);
b_out  = repmat(beta_hat0, [1 1 K]);
a_out  = zeros(K, p);

b_old = ones(N, p);
a_old = zeros(1, p);

cvx_quiet(true)

cvx_solver mosek 
% this line can be commented out if MOSEK is not installed
% the CVX default solver runs much more slowly than mosek

for r = 1:R
    for k = 1:K
        
        
        for kk = setdiff(1:K, k)
            pen(:, kk) =  norms( bsxfun(@minus, b_out(:, :, kk), a_out(kk, :) ), 2, 2 );
        end
        pen(:,k) = ones(N, 1);
        penalty_out  =  prod(pen, 2);
        
        cvx_begin
        variable b(N, p)
        variable a(1, p);
        B1 = kron( b, ones(T,1) );
        Q =  1/(N*T) * sum_square(  y - sum(X .* B1, 2) );
        
        
        pen_k = norms(  b - repmat(a, N, 1), 2, 2 );
        penalty =  penalty_out' * pen_k ;
        
        
        minimize( Q  + lambda/N *  penalty);
        cvx_end
        pen(:,k) = pen_k;
        
        b_out(:, :, k) = b;
        a_out(k, :) = a;
        
    end
    
    b_new = b_out(:, :, k);
    a_new = a_out(k, :);
    
    
    if criterion( a_old, a_new, b_old, b_new, tol  ) == 1
        break;
    end
    
    
    a_old = a_out(k, :);
    b_old = b_out(:, :, k);
    
end

% put b_out into the nearest a_out

    distance = zeros(N,K);
    group_est = zeros(N, K);
    for k = 1:K

        nominator = norms( bsxfun(@minus, b_out(:, :, k), a_out(k, :) ), 2 ,2 );
        distance(:,k) = nominator./norm( a_out(k,:) ) ;

    end
    [~, which] = min(distance, [], 2);
    b_est = 999 * ones(N, p );

    for i = 1:N
        group_est(i, which(i)) = 1;
        b_est(i, :)  = a_out(which(i), :);
    end

end

function [d] = criterion( a_old, a_new, b_old, b_new, tol)

    d = 0;
    a_nominator = sum( abs( a_old - a_new ) );
    a_denominator = sum( abs( a_old ) ) + 0.001;

    b_nominator = mean( abs( b_old - b_new ) );
    b_denominator = mean( abs( b_old ) ) + 0.001;

    if a_nominator/ a_denominator < tol &&  b_nominator / b_denominator < tol
        d = 1;
    end

end