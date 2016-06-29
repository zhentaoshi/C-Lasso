function [ group_binary ] = groups_coerce(b, a, K)



        N = size(b, 1);
        distance = zeros(N, K);
        for k = 1:K
            distance(:,k) = norms( bsxfun(@minus, b, a(k,:) ), 2, 2);
        end
        [ ~, group] = min(distance,[], 2);

        group_binary = zeros( size(group) );
        for k = 1:K
            group_binary(:,k) = (group == k);
        end
end
