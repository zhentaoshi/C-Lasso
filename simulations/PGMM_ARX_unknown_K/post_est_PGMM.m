function [a_corr, post_a, post_a_corr, vari_a, vari_a_post] = post_est_PGMM( N, T, a_hat, dat, this_group )

global p d
D = p + d - 1;
    index = 1:N;
    g_index = index(this_group);
    g_data = dat( ismember(dat.N, g_index), : );

    y = g_data.y;
    X = g_data.X;
    Z = g_data.Z;

    if all( logical( this_group ) == 0 )
        post_a = zeros(p, 1);
        a_corr = post_a;
        post_a_corr = post_a;
		    vari_a = ones(1,p);
    vari_a_post = ones(1,p);
	elseif sum(this_group) >= p

        XZZX = (X' * Z) * (Z' * X);
        XZZy = (X' * Z) * (Z' * y);
        post_a = (pinv(XZZX) * XZZy)';

        a_corr = a_hat' - bias_PGMM_ARX(T, a_hat', y, X, Z );
        post_a_corr = post_a' - bias_PGMM_ARX(T, post_a', y, X, Z );

		    [vari_a] = var_post_GMM(T, a_hat', y, X, Z);
    [vari_a_post] = var_PGMM_ARX(T, post_a', y, X, Z);
    vari_a = diag(vari_a)';
    vari_a_post = diag(vari_a_post)';

    else

        a_corr = a_hat;
        XZZX = (X' * Z) * (Z' * X);
        XZZy = (X' * Z) * (Z' * y);
        post_a = (pinv(XZZX) * XZZy)';
        post_a_corr = post_a;

		    [vari_a] = var_post_GMM(T, a_hat', y, X, Z);
    [vari_a_post] = var_PGMM_ARX(T, post_a', y, X, Z);
    vari_a = ones(1,p);
    vari_a_post = ones(1,p);
    end
end
