function [bia] = bias_PGMM(T, a, y, X, Z)














global p d
D = p+d-1;

MT = ceil(2 * T^(1/4));
Nk = size(y, 1)/T;



Q =  (1/T) * Z' * X;
A = 1/(Nk) * (Q' * Q);

u = y - X * a;


uu = zeros(Nk, T);
XX = zeros(Nk, p, T);
ZZ = zeros(Nk, D, T);

for i = 1:Nk
    uu( i, :) = u( ( T*(i-1) + 1 ): (T*i) );
    XX( i, :, :) = X( ( T*(i-1) + 1 ): (T*i), : )';
    ZZ( i, :, :) = Z( ( T*(i-1) + 1 ): (T*i), : )';
end

kernel = @( u ) (1- u / MT ) .* ( ( u <= MT ) * (MT >=0 ) ) / ((MT+1)/2) ;

combNK = combnk(1:T, 2);
comb_st = [ repmat( (1:T)', 1, 2  ); combNK; combNK( :, [2,1] )  ];

TTT = size(comb_st, 1);
kern_st = zeros(1, TTT);
XZZU_st = zeros(p, TTT);
for j = 1:TTT
    s = comb_st(j,1);
    t = comb_st(j,2);
    kern_st(j) = kernel( s - t );

     XXis = XX( :, :, s);
     ZZis = ZZ( :, :, s);
     ZZit = ZZ( :, :, t);
     uuit = uu( :, t);

     XZZU_st(:, j) = (XXis' * ZZis * ZZit' * uuit)';
end

kXZZU = sum( bsxfun(@times, XZZU_st, kern_st), 2);
B = 1/(Nk^0.5 * T^(1.5) ) * kXZZU;
bia =  -1/(Nk^0.5 * T^0.5 ) * pinv(A) * B;
end
