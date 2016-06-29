function [class] = report(hat_group, hat_group_coerce, K)



global N_cut

N = size(hat_group, 1);

group0 = zeros(N, K);
classified = logical( sum(hat_group, 2) );


for k = 1:K
    if k == 1
        group0(1: N_cut(k) , k) = 1;
    else
        group0( (N_cut(k-1) + 1):N_cut(k), k ) = 1;
    end
end
group0 = logical(group0);
proportion = [0.3 0.3 0.4];

hat_group_classified = hat_group(classified, :);
group0_classified    = group0(   classified, :);


Eik = zeros(1, K);
Fik = zeros(1, K);
for k = 1:K
    Eik(k) = mean( ~hat_group_classified( group0_classified(:,k  ),k ) ) ;
    Fik(k) = mean( hat_group_classified(~group0_classified(:,k ), k ) );
end
Ei = Eik * proportion';
Fi = Fik * proportion';

miss = 1-mean( sum(hat_group .*group0, 2) );




hat_group_coerce = logical(hat_group_coerce);

Eik_co = zeros(1,K);
Fik_co = zeros(1,K);
for k = 1:K
    Eik_co(k) =  mean (~hat_group_coerce(group0(:,k), k)  );
    Fik_co(k) =  mean( group0( ~hat_group_coerce(:,k), k ) );
end
Ei_co = Eik_co * proportion';
Fi_co = Fik_co * proportion';

miss_co = 1-mean( sum(hat_group_coerce .* group0, 2) );


class = [Ei, Fi, miss, Ei_co, Fi_co, miss_co];
class
end
