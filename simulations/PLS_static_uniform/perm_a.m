function [perm_hat ] = perm_a(a0, a)




    K = size(a0, 1);
    order_perm = perms(1:K);
    length_perm = size(order_perm, 1);

    dist = zeros(1, length_perm);
    for i = 1:length_perm
        a_perm= a( order_perm(i, :), :);

        dist(i) = norm(a0-a_perm, 2);
    end
    [~, which] = min(dist);
    perm_hat = order_perm(which, :);

end
