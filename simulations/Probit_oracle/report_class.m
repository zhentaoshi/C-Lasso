function [underclass, correctclass, wrongclass] = ...
    report_class( b, K, est_group )















    global a0 N_cut

    N = size(est_group,1);
    index_n = 1:N;
    class_index = index_n( sum(est_group, 2) == 1 );
    under_class = mean( sum(est_group, 2) ~= 1 );


    class = 999 * ones(1,N);
    for k = 1:K
        if k == 1
            range = 1:N_cut(1);
        else
            range = (N_cut(k-1) + 1): N_cut(k);
        end


        for i = intersect( range, class_index )
            bi = reshape( b(i,:), 2, 3);
            dist = norms( bi' - a0 , 2 , 2);
            if all( dist(k) <= dist )
                class(i) = 1;
            else
                class(i) = 0;
            end
        end
    end

    underclass = under_class;
    correctclass = mean( class == 1);
    wrongclass   = mean( class == 0);



end
