function [vari] = var_post_GMM(T, a, y, X, Z)















d = size(Z,2);
p = size(X,2);
MT = floor(  T^(1/4));
Nk = size(y, 1)/T;


u = y - X * a ;

W = pinv( Z' * Z/(Nk * T) );
Q = X' * Z/(Nk * T);
A = Q * W * Q';


kernel = @( u ) (1- u / MT ) .* ( ( u <= MT ) * (MT >=0 ) )  ;

combNK = combnk(1:T, 2);
comb_st = [ repmat( (1:T)', 1, 2  ); combNK; combNK( :, [2,1] )  ];

TTT = size(comb_st, 1);
kern_st = zeros(1, 1, TTT);
ZuuZ_st = zeros(d,d, TTT);

ZZ = reshape(Z, [T Nk d] );
ZZ = permute(ZZ, [2 3 1]);
uu  = reshape( u, [T Nk])';
for j = 1:TTT
    s = comb_st(j,1);
    t = comb_st(j,2);
    kern_st(j) = kernel( s - t );

     Zus =   ZZ(:, :, s) .* repmat(uu(:, s), [1 d]) ;
     Zut =   ZZ(:, :, t) .* repmat(uu(:, t), [1 d]) ;

     Zus = repmat( Zus, [1 1 d]);

     Zut = permute( Zut, [1 3 2]);
     Zut = repmat( Zut, [1 d 1]);

     ZuuZ_st(:, :, j) = permute( sum( Zus.* Zut, 1), [2 3 1] );
end

kZuuZ = sum( ZuuZ_st .* repmat(kern_st, [d d 1] ), 3)/(Nk*T);
V = Q * W * kZuuZ * W * Q';

vari =  pinv(A) * V * pinv(A)/(Nk * (T-1) );
end
