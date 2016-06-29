function [a] = post_est( N, T, dat, this_group, bias_corr )

global p
    this_group = logical(this_group);
    index = 1:N;
    g_index = index(this_group);
    g_data = dat( ismember(dat.N, g_index), : ); 
    
    if isempty(this_group)
        a = zeros(p, 1);
    elseif sum(this_group) >= p
        a = regress( g_data.y, g_data.X );
		if bias_corr == 1
			a = a - bias_PLS( T, a, g_data.y, g_data.X );
		end
    else 
        XX = g_data.X' * g_data.X;
        Xy = g_data.X' * g_data.y;
        a = pinv(XX) * Xy;
    end
   
end