function [vari] = var_PGMM_ARX(T, a, y, X, Z)















D = size(Z,2);
p = size(X,2);
MT = ceil(2 * T^(1/4));
Nk = size(y, 1)/T;


u = ( y - X * a );


uu = zeros(Nk, T);
XX = zeros(Nk, p, T);
ZZ = zeros(Nk, D, T);
QZ = zeros(Nk, p, D);
Q_bar = zeros(p, p, Nk);

for i = 1:Nk
    uu( i, :) =  u( ( T*(i-1) + 1 ): (T*i) ) ;
    XX( i, :, :) = X( ( T*(i-1) + 1 ): (T*i), : )';
    ZZ( i, :, :) = Z( ( T*(i-1) + 1 ): (T*i), : )';
    XXi = X( ( T*(i-1) + 1 ): (T*i), : )';
    ZZi = Z( ( T*(i-1) + 1 ): (T*i), : )';
    QZ( i, :, :) = XXi * ZZi' / T;
    Q_bar(:, :, i) = (XXi * ZZi') * (ZZi * XXi') / T^2;
end
A_bar = mean(Q_bar, 3);


kernel = @( u ) (1- abs(u) / MT ) .* (abs(u) <= MT);

combNK = combnk(1:T, 2);
comb_st = [ repmat( (1:T)', 1, 2  ); combNK; combNK( :, [2,1] )  ];

TTT = size(comb_st, 1);
kern_st = zeros(1, TTT);
XZZU_st = zeros(p, TTT);
QZUUZQ_st = zeros(p,p,TTT);

for j = 1:TTT
    s = comb_st(j,1);
    t = comb_st(j,2);
    kern_st(j) = kernel( s - t );

     XXis = XX( :, :, s);
     ZZis = ZZ( :, :, s);
     ZZit = ZZ( :, :, t);
     uuit = uu( :, t) ;
     uuis = uu( :, s) ;

     XZZU_st(:, j) = (XXis' * ZZis * ZZit' * uuit)';


     Zuit =  ZZit.*repmat( uuit, [ 1 D]);
     Zuit = permute(Zuit, [1 3 2]);
     Zuit = repmat(Zuit, [1 p 1]);
     QZuit  = sum( QZ .* Zuit, 3);
     QZuit = repmat( QZuit, [1 1 p]);

     Zuis = ZZis .* repmat( uuis, [1 D]);
     Zuis = permute( Zuis, [1 3 2]);
     Zuis = repmat( Zuis, [1 p 1]);
    QZuis  = sum( QZ .* Zuis, 3);
    QZuis = permute(QZuis, [1 3 2]);
    QZuis = repmat( QZuis, [1 p 1]);

    QZUUZQ_st(:,:, j)  = permute( sum( QZuis .* QZuit , 1 ), [ 2 3 1]);
end
kern_st = kern_st/ max(kern_st);
kXZZU = sum( bsxfun(@times, XZZU_st, kern_st), 2);
B = 1/(Nk^0.5 * T^(1.5) ) * kXZZU;


kern_mat = permute(kern_st, [3 1 2]);
kern_mat = repmat(kern_mat, [p p 1]);
kQZUUZQ = sum( QZUUZQ_st.*kern_mat, 3);
V =  kQZUUZQ/(Nk * T ) ;

vari =  pinv(A_bar) * V * pinv(A_bar)/(Nk * T);
end
