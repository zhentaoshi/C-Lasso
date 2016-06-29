function [a_corr, post_a, post_a_corr, vari_a, vari_a_post] = post_est_PLS( N, T, a_hat, dat, this_group, bias_corr )

global p
index = 1:N;
g_index = index(this_group);
g_data = dat( ismember(dat.N, g_index), : );

if all( logical( this_group ) == 0 )
    post_a = zeros(p, 1);
    a_corr = post_a;
    post_a_corr = post_a;

elseif sum(this_group) < p

    post_a = a_hat;
    a_corr = a_hat;
    post_a_corr = a_hat;
else

    post_a = regress( g_data.y, g_data.X );
    if bias_corr == 1
        a_corr = a_hat' - bias_PLS(T, a_hat', g_data.y, g_data.X );
        post_a_corr = post_a - bias_PLS(T, post_a, g_data.y, g_data.X );
    else
        a_corr = a_hat;

        XX = g_data.X' * g_data.X;
        Xy = g_data.X' * g_data.y;
        post_a = pinv(XX) * Xy;
        post_a_corr = post_a;
    end
end


[vari_a] = var_PLS(T, a_hat', g_data.y, g_data.X );
[vari_a_post] = var_PLS(T, post_a, g_data.y, g_data.X );
vari_a = diag(vari_a)';
vari_a_post = diag(vari_a_post)';

end
