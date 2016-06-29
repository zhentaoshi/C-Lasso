function [b_out, a_out, c_out] = PNL_est(N, T, b_initial,  y, X, K, lam, R, tol )
% estimation. Iterate the bi-convex function

% INPUT:
%   y, X1, X2: the data
%   R: maximum number of iterations
%   tol: tolerance level

global p

pen  = ones(N, K);
b_out  = repmat( b_initial, [ 1 1 K]) ;
a_out  = zeros(K, p);

a_old = zeros(K, p);
optval_old = 999;
sign_y = (2*y - 1);

cvx_quiet(true)
%%

if K == 1 % no classification is needed.
    cvx_begin
    %    cvx_solver mosek
    variable c(N, 1);
    variable b(1, p);
    
    cst = kron(c, ones(T,1) );
    Q =  -1/(N*T) * sum( log_normcdf( ( cst + X * b' ) .* sign_y ) );
    
    minimize(Q);
    cvx_end
    
    a_out = b;
    b_out = b;
    c_out = c;
else
    for r = 1:R
        QQQ = 999*ones(K,1);
        for k = 1:K
            
            
            % calculate the fixed part of the penalty
            for kk = setdiff(1:K, k)
                pen(:, kk) =  norms( bsxfun(@minus, b_out(:, :, kk), a_out(kk, :) ), 2, 2 );
            end
            pen(:,k) = ones(N, 1);
            penalty_out  =  prod(pen, 2);
            
            cvx_begin
            % cvx_solver mosek
            
            variable c(N, 1);
            variable b(N, p);
            variable a(1, p);
            
            B1 = kron( b, ones(T,1) );
            cst = kron(c, ones(T,1) );
            Q =  -1/(N*T) * sum( log_normcdf( ( cst + sum(X .* B1, 2) ) .* sign_y ) ); % + lam * norm(c,1);
            %
            %  penalty
            pen_k = norms(  b - repmat(a, N, 1), 2, 2 );
            penalty =  penalty_out' * pen_k ;
            %
            %  objective
            QQ = Q  + lam/N *  penalty;
            minimize( QQ);
            cvx_end
            QQQ(kk) = QQ;
            
            pen(:,k) = pen_k;
            
            b_out(:, :, k) = b;
            a_out(k, :) = a;
            c_out = c;
        end % the K iteration
        % at the end of the iteration
        % update the parameter
        % test the convergence criterion
        
        a_new = a_out;
        optval_new = sum(QQQ);
        
        if criterion2( a_old, a_new, optval_old, optval_new, tol  ) == 1
            break;
        end
        
        % update the parameter
        a_old = a_new;
        optval_old = optval_new;
        
    end
end
end