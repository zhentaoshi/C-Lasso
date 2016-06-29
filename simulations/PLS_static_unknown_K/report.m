function [class] = report( hatgroup, K)



global N_cut

N = size(hatgroup, 1);

group0 = zeros(N, K);


for k = 1:K
    if k == 1
        group0(1: N_cut(k) , k) = 1;
    else
        group0( (N_cut(k-1) + 1):N_cut(k), k ) = 1;
    end
end
group0 = logical(group0);
hatgroup = logical(hatgroup);

Eik = zeros(1,K);
Fik = zeros(1,K);
for k = 1:K
    Eik(k) =  mean (~hatgroup(group0(:,k), k)  );
    Fik(k) =  mean( group0( ~hatgroup(:,k), k ) );
end
Ei_co = Eik_co * proportion';
Fi_co = Fik_co * proportion';


class = [0,0,0, Ei_co, Fi_co, 0];
class;
end
