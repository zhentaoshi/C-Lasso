function [b_out, a_out] = PLS_est(N, T, y, X, b0, K, lambda, R, tol)
    % Su, Shi and Phillips (2016)
    % estimate the penalized least square.
    % this is the core function that performs the optimization
    
    % INPUT:     
    %   N: sample size
    %   T: number of time periods
    %   y: dependent varaible
    %   X: independent variable(s)
    %   b0: a vector of initial estimates. size(b0) = N * # of coefficients
    %   K: number of groups. "K = 1" is allowed
    %   R: maximum number of iterations
    %   tol: tolerance level
    
    % OUTPUT:
    %   b_out: beta estimate
    %   a_out: alpha estimate
   
    p = size(X, 2);
    
    pen  = ones(N, K);
    b_out  = repmat(b0, [ 1 1 K]);
    a_out  = zeros(K, p);
    
    a_old = 999 * zeros(K, p);
    optval_old = 999;
    
    cvx_quiet(true)
    %%
    
    if K == 1
        cvx_begin
            variable b(N, p)
            variable a(1, p);

            B1 = kron( b, ones(T,1) );
            Q =  1/(N*T) * sum_square(  y - sum(X .* B1, 2) );

            % the penalty
            penalty = sum( norms(  b - repmat(a, N, 1), 2, 2 ) );

            % objective
            minimize( Q  + lambda/N *  penalty ) ;
        cvx_end
        
        a_out = a;
        b_out = b;
    else
        for r = 1:R
            
            QQQ = 999 * ones(1,K);
            for k = 1:K            

                % calculate the fixed part of the penalty
                for kk = setdiff(1:K, k)
                    pen(:, kk) =  norms( bsxfun(@minus, b_out(:, :, kk), a_out(kk, :) ), 2, 2 );
                end
                pen(:,k) = ones(N, 1);
                penalty_out  =  prod(pen, 2);

                cvx_begin
                cvx_solver sedumi
                    variable b(N, p)
                    variable a(1, p);
                    B1 = kron( b, ones(T,1) );
                    Q =  1/(N*T) * sum_square(  y - sum(X .* B1, 2) );

                    % the penalty
                    pen_k = norms(  b - repmat(a, N, 1), 2, 2 ) ;
                    penalty =  penalty_out' * pen_k;

                    % objective
                    QQ = Q  + lambda/N *  penalty;
                    minimize( QQ );    
                    subject to
                        a(1) <= 0.95; 
                        a(1) >= 0;
                cvx_end
                pen(:,k) = pen_k; % update the penalty

                b_out(:, :, k) = b;
                a_out(k, :) = a;
                
                QQQ(k) = QQ;
            end % the K iteration
            % at the end of the iteration
            % update the parameter
            % test the convergence criterion
            optval_new = sum(QQQ);
            [r optval_new, penalty];
            a_new = a_out;

            if criterion2( a_old, a_new, optval_old, optval_new, tol  ) == 1 
                break;
            end

            % update the parameter
            a_old = a_new;
            optval_old = optval_new;

        end
    end
end