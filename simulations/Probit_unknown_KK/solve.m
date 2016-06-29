function [b_out, c_out] = solve(a_initial, Nk, T, X_vector, y_vector)


probit_in_b = @(b) probit( b, Nk, T, X_vector, y_vector);

tolf = 1e-8;
optionscon = optimset('Algorithm', 'trust-region-reflective',...
    'Display', 'off', ...
    'TolFun', tolf, 'TolX', tolf, 'GradObj', 'on', 'Hessian', 'on',...
    'MaxIter', 1e+2);

	sign_y_all = 2*y_vector - 1;


probit_in_b = @(b) probit( b, Nk, T, X_vector, y_vector);



[b_out, fval] = fminunc(probit_in_b, a_initial,optionscon);

    cvx_begin
        variable c(Nk, 1);
        cst = kron(c, ones(T,1) );
        Q =  - 1/(Nk*T) * sum( log_normcdf( ( cst + X_vector * b_out ) .* sign_y_all ) );
        minimize( Q );
    cvx_end

	c_out = c;
end


