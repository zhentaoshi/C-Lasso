function [bia] = bias_PLS(T, a, y, X )


























global p
MT = ceil( T^(1/6));
Nk = size(y, 1)/T;



Phi = 1/(Nk * T) * (X' * X);


Xa = X * a;
ui = y - Xa;


uis = zeros( T, Nk);
XX = zeros( T, p, Nk);

for i = 1:Nk
    u_within_group = ui( ( T*(i-1) + 1 ): (T*i) );
    uis( :, i) = bsxfun(@minus, u_within_group, mean(u_within_group));
    XX(:, :, i) = X( ( T*(i-1) + 1 ): (T*i), : );
end


kernel = @( u ) (1- u / MT ) .* ( ( u <= MT ) * (u >=0))  /(( MT+1)/2) ;


combNK = combnk(1:T, 2);
comb_st = [ repmat( (1:T)', 1, 2  ); combNK; combNK( :, [2,1] )  ];

TTT = size(comb_st, 1);
kern_st = zeros(1, TTT);
XU_st = zeros(p, TTT);

for j = 1:TTT
    s = comb_st(j,1);
    t = comb_st(j,2);
    kern_st(j) = kernel( s-t );

    XXis = XX( s, :, :);
    XXis = permute(XXis, [2 3 1]);
    uuit = uis( t, :);

    XU_st(:, j) =   XXis * uuit';
end

kXU = sum( bsxfun(@times, XU_st, kern_st), 2);
B = -1/(Nk^0.5 * T^(1.5) ) * kXU;
bia = (1/(Nk * T)^0.5) * pinv(Phi) * B;

end
