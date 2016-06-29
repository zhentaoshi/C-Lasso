function sd_rho = SD_rho_robust( hat_a, fe, N, T, y, X)


p = 2;
MT = ceil( T^(1/4));


YL = y; YR = X(:,:,1); X = X(:,:,2);

[logl, grad, Hess,  Grho, Gbeta]=G_H(hat_a, fe,YL,YR,X);

Grho_vector = reshape(Grho, [T*N 1]);
Gbeta_vector= reshape(Gbeta,[T*N 1]);
Nk = size(YL, 2);




Grad = [Grho_vector, Gbeta_vector];

GG = zeros( T, p, Nk);
for i = 1:Nk
    G_within_group = Grad( ( T*(i-1) + 1 ): (T*i), : );
    GG( :, :, i) = G_within_group;

end


kernel = @( u )  ( abs(u) <= MT );
combNK = combnk(1:T, 2);
comb_st = [ repmat( (1:T)', 1, 2  ); combNK; combNK( :, [2,1] )  ];

TTT = size(comb_st, 1);
kern_st = zeros(1, 1, TTT);
GG_st = zeros(p, p, TTT);


for j = 1:TTT
    s = comb_st(j,1);
    t = comb_st(j,2);
    kern_st(j) = kernel( s-t );

    Gis = GG( s, :, :);
    Gis = permute(Gis, [3 2 1]);

    Git = GG( t, :, :);
    Git = permute(Git, [3 2 1]);

    GG_st(:,:,j) = Gis' * Git;
end
kern_mat = repmat( kern_st, [p p 1]);
Omega = sum( GG_st .* kern_mat, 3 )/(T * Nk);


vari = pinv(-Hess) * Omega * pinv(-Hess);
sd_rho = sqrt( vari(1,1)/(N*(T-1) - 2 ) );


end



