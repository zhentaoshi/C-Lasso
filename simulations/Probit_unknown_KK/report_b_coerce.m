function [a_out, beta_est, group] = report_b_coerce( b, a, K )










global a0 p
a_true = a0;





K_true = size(a0,2);
if K > K_true;
    a_true = [ a_true; zeros(K - K_true, p) ];
end


if K == 1
    beta_est = b;
    group = 1;
else


dist_KK = zeros(K,K);
for k1 = 1:K
    for k2 = 1:K
        dist_KK(k1,k2) = sum( (a(k1,:) - a_true(k2,:) ).^2, 2);
    end
end

order_M = zeros(K, 2);
order_M(:,1) = 1:K;

[~, choice1] = min(min(dist_KK));
[~, location1] = min(dist_KK(:, choice1));
order_M(choice1,2) = location1;

dist_KK(:, choice1) = 9999;
dist_KK(location1, :) = 9999;
[~, choice2] = min(min(dist_KK));
[~, order_M(choice2,2)] = min(dist_KK(:, choice2));

order_M( order_M(:,2) == 0, 2) =  setdiff(1:K, order_M(:,2) );




    a = a(order_M(:, 2), :);

    N = size(b, 1);
    p = size(b,2);

    group = zeros(N, 1);
    beta_est = zeros(N, p);
    for i = 1:N
        b_i = b(i, :, :);
        b_i = permute(b_i, [3 2 1]);
        dist_i = sum( (b_i - a).^2, 2 );
        [~, group(i)] = min(dist_i);
        beta_est(i,:) = a( group(i),:);
        a_out = a;
    end

end
end
