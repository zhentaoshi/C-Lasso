function [b_out, a_out] = PGMM_est(y, X, Z, b0, W, K, lambda, R, tol)
    % The corrct moment conditions is 
    % E[ (y - Xb) * Z ] = 0
    
    % INPUT: 
    %   y: dim = T * N
    %   X: dim = T * N * p
    %   Z: dim = T * N * (d1 + d2 + ... + dp). where d1 is the dimension of
    %   the instruments for p = 1, d2 is for p = 2, and so on.
    %   W: dim = (Nd) * (Nd). Each diagoal block is d*d.
    
    %   R: maximum number of iterations
    %   tol: tolerance level
    
    % OUTPUT:
    
    global N T p d 
    
    pen  = ones(N, K);
    b_out = repmat(b0, [1 1 K]);
    
    
    a_out  = zeros(K, p);
    a_old = zeros(K, p);
    
    optval_old = 999;
    cvx_quiet(true)
    %%
    
    if K == 1
        a_out = regress(y, X);
        a_out = a_out';
        b_out(:,:,K) = repmat( a_out, N, 1);
    else
        for r = 1:R
            for k = 1:K            

                % calculate the fixed part of the penalty
                for kk = setdiff(1:K, k)
                    pen(:, kk) =  norms( bsxfun(@minus, b_out(:, :, kk), a_out(kk, :) ), 2, 2 );
                end
               pen(:,k) = 1; 
                penalty_out  =  prod(pen, 2);

                cvx_begin
                cvx_solver mosek
                    variable b(N, p)
                    variable a(1, p);

                    B1 = repmat( permute( b, [3 1 2] ), [T 1 1] );
                    e = y - sum( X .* B1, 3); % sum over the p dim. size(e) = T * N

                    g = repmat(e, [1 1 size(Z,3) ] ) .* Z ; % size(g) = T * N * d
                    g = sum(g)/size(g,1); 
                    g = reshape( g, N*size(Z,3) , 1 ); %make sure the dimension are right.
                    Q = (1/N) * quad_form(g', W ); % each diagonal block is d*d.

                    % the penalty
                    pen_k = norms(  b - repmat(a, N, 1), 2, 2 );
                    penalty =  penalty_out' * pen_k ;

                    % objective
                    QQ = Q  + lambda/ N *  penalty;
                    minimize( QQ );    
                     subject to 
                       % a(1) <= .95; % this restriction is particularly impose here
                       % a(1) >= 0;
                cvx_end
                pen(:,k) = pen_k; % update the penalty

                b_out(:, :, k) = b;
                a_out(k, :) = a;
                QQQ(k) = QQ;
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

        end % the R iteration
    end
    
    
    
    
end