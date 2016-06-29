function [bia] = bias_PGMM_ARX(T, a, y, X, Z)
% this function calculate A and B, the two components of the bias term as
% in the note p.2.
%
% INPUT
% y, X, Z are not those from the full dataset. They are from one group.
% a is the post-estimator without bias correction.
%
% OUTPUT
%
% in the function "post_est_PGMM", 
% size(Z) = (Nk * T) * d;
% size(X) = (Nk * T) * p;
% size(y) = (Nk * T) * 1;
% demean = @(s) s - mean(s);

d = size(Z,2);
p = size(X,2);
MT = ceil(2 * T^(1/4)); 
Nk = size(y, 1)/T; 


%%
u = ( y - X * a );
uu = zeros(Nk, T);
XX = zeros(Nk, p, T);
ZZ = zeros(Nk, d, T);

Qi = zeros( p, p);
for i = 1:Nk
    uu( i, :) =  u( ( T*(i-1) + 1 ): (T*i) ) ;
    Xi = X( ( T*(i-1) + 1 ): (T*i), : )'; %  p * T
    XX( i, :, :) = Xi;
    
    
    Zi= Z( ( T*(i-1) + 1 ): (T*i), : )'; % d * T
    ZZ( i, :, :)  = Zi; 
    Qi = Qi +  (Xi * Zi' * Zi * Xi'/T^2)   ; % p * p
end
A_bar = Qi / Nk;
%%
kernel = @( u ) (1- abs(u) / MT ) .* (abs(u) <= MT);

combNK = combnk(1:T, 2);
comb_st = [ repmat( (1:T)', 1, 2  ); combNK; combNK( :, [2,1] )  ]; % all combinations, including self

TTT = size(comb_st, 1);
kern_st = zeros(1, TTT);
XZZU_st = zeros(p, TTT);

for j = 1:TTT
    s = comb_st(j,1);
    t = comb_st(j,2);
    kern_st(j) = kernel( s - t );
    
     XXis = XX( :, :, s);
     ZZis = ZZ( :, :, s);
     ZZis = permute( ZZis, [1 3 2]);     
     XZis = repmat( XXis, [1 1 d]) .* repmat(ZZis, [1 p 1]);

     ZZit = ZZ( :, :, t);
     uuit = ( uu( :, t) );
     Zuit = ZZit .* repmat( uuit, [1 d]);
     Zuit = permute(Zuit, [1 3 2]);
     Zuit = repmat( Zuit, [1 p 1]);
     
     XZZu1 = sum( XZis .* Zuit, 3);
     XZZu1 = sum( XZZu1, 1);
     XZZU_st(:, j) = XZZu1';
end

kXZZU = sum( bsxfun(@times, XZZU_st, kern_st), 2); 
B = 1/(Nk^0.5 * T^(1.5) ) * kXZZU;

bia =  1/(Nk^0.5 * T^0.5 ) * pinv(A_bar) * B;
end