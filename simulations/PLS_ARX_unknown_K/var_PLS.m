function [vari] = var_PLS(T, a, y, X )












global p
MT = ceil( 2* T^(1/4));
Nk = size(y, 1)/T;


Phi = 1/(Nk * T) * (X' * X);

Xa = X * a;
ui = y - Xa;

uis = zeros( T, Nk);
XX = zeros( T, p, Nk);

for i = 1:Nk
    u_within_group = ui( ( T*(i-1) + 1 ): (T*i) );
    uis( :, i) = u_within_group;
    XX(:, :, i) = X( ( T*(i-1) + 1 ): (T*i), : );
end


kernel = @( u ) (1- abs(u) / MT ) .*  ( abs(u) <= MT );
combNK = combnk(1:T, 2);
comb_st = [ repmat( (1:T)', 1, 2  ); combNK; combNK( :, [2,1] )  ];

TTT = size(comb_st, 1);
kern_st = zeros(1, 1, TTT);
XuuX_st = zeros(p, p, TTT);

for j = 1:TTT
    s = comb_st(j,1);
    t = comb_st(j,2);
    kern_st(j) = kernel( s-t );

    XXis = XX( s, :, :);
    XXis = permute(XXis, [3 2 1]);
    uuis = uis( s, :);
    uuis = repmat( uuis', [1 p]);
    Xuis = XXis .* uuis;
    Xuis = repmat( Xuis, [1 1 p]);

    XXit = XX( t, :, :);
    XXit = permute(XXit, [3 2 1]);
    uuit = uis( s, :);
    uuit = repmat( uuit', [1 p]);
    Xuit = XXit .* uuit;
    Xuit = permute( Xuit, [1 3 2]);
    Xuit = repmat( Xuit, [1 p 1]);

    XuuX_st1 = sum( Xuis .* Xuit, 1);
    XuuX_st(:,:, j) =  permute( XuuX_st1, [2 3 1]);
end
kern_mat = repmat( kern_st, [p p 1]);
kXuuX = sum( XuuX_st .* kern_mat, 3 )/(T * Nk);


vari = 1/(Nk * T) * pinv(Phi) * kXuuX * pinv(Phi);
end
