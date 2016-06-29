function [b_out, a_out, opt_r] = PLS_est_Q(N, T, y, X, br, b0, K, lambda, R, tol)
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
b_out  = repmat(br, [ 1 1 K]);
b_out(:,:,1) = b0;
a_out  = zeros(K, p);
    
b_old = ones(N, p);
a_old = zeros(K, p);

    optval_old = 999;
    
    opt_r = zeros(80,1);
    cvx_quiet(true)
    %%
    

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
                    penalty =  penalty_out' * pen_k ;

                    % objective
                    QQ = Q  + lambda/N *  penalty;
                    minimize( QQ );    
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
            opt_r(r) = optval_new;
            a_new = a_out;

            if criterion2( a_old, a_new, optval_old, optval_new, tol  ) == 1 
                break;
            end

            % update the parameter
            a_old = a_new;
            optval_old = optval_new;

        end
end