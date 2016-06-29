function [a_corr, post_a, post_a_corr] = post_est_PGMM( N, T, a_hat, dat, this_group )
% all data is restricted to that group.

    index = 1:N;
    g_index = index(this_group);
    g_data = dat( ismember(dat.N, g_index), : ); % group-specific data
    
    y = g_data.y; %note the change of definition. This is correct.
    X = g_data.X;
    Z = g_data.Z;
    

    

    


        XZZX = (X' * Z) * inv(Z' * Z) * (Z' * X);
        XZZy = (X' * Z) * inv(Z' * Z) * (Z' * y);
        post_a = (pinv(XZZX) * XZZy)';

            a_corr = a_hat' - bias_PGMM_ARX(T, a_hat', y, X, Z);
			post_a_corr = post_a' - bias_PGMM_ARX(T, post_a', y, X, Z );


end