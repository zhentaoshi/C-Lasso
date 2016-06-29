function theta_bar = SPJ_zt(a_initial, y_vector, X_vector, Nk, T)

global p

period1_i = [ones(floor(T/2), 1); zeros( ceil(T/2),1)];
period1 = logical( kron(  ones(Nk, 1), period1_i )) ;
period2 = logical( 1-period1);

X_vector1 = X_vector(period1, :);
y_vector1 = y_vector(period1);

X_vector2= X_vector(period2, :);
y_vector2 = y_vector(period2);


a_tt = zeros(1,2);
for tt = 1:2

    if tt == 1
	TT = sum(period1_i);
    else
	TT = T - sum(period1_i);
	end

    if tt == 1
        X_vector_pp = X_vector1;
        y_vector_pp = y_vector1;
    elseif tt == 2
        X_vector_pp = X_vector2;
        y_vector_pp = y_vector2;
    end

    [b] = solve(a_initial, Nk, TT, X_vector_pp, y_vector_pp);

    a_tt(tt) = b(1);
cvx_begin
end
theta_bar =  mean(a_tt);
end
