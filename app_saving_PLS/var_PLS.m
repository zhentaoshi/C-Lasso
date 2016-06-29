function [vari] = var_PLS(T, a, y, X )
% Liangjun Su, Zhentao Shi, Peter Phillips (2014)

%%
% for LS, the variance of C-lasso and post Lasso are the same. 
%
% INPUT
%   y, X, Z are not those from the full dataset. They are from one group.
%   a: the post-estimator without bias correction.
%
% OUTPUT
%   vari: the variance of the bias-corrected estimator

global p
MT = ceil( T^(1/4));
Nk = size(y, 1)/T; 

%% 
Phi = 1/(Nk * T) * (X' * X); % correct

Xa = X * a;
ui = y - Xa; 
%% 
uis = zeros( T, Nk);
XX = zeros( T, p, Nk);

for i = 1:Nk
    u_within_group = ui( ( T*(i-1) + 1 ): (T*i) );
    uis( :, i) = u_within_group; % bsxfun(@minus, u_within_group, mean(u_within_group));
    XX(:, :, i) = X( ( T*(i-1) + 1 ): (T*i), : );   
end

%%
kernel = @( u ) (1- abs(u) / MT ) .*  ( abs(u) <= MT );
combNK = combnk(1:T, 2);
comb_st = [ repmat( (1:T)', 1, 2  ); combNK; combNK( :, [2,1] )  ]; % all combinations, including self

TTT = size(comb_st, 1);
kern_st = zeros(1, 1, TTT);
XuuX_st = zeros(p, p, TTT);

for j = 1:TTT
    s = comb_st(j,1);
    t = comb_st(j,2);
    kern_st(j) = kernel( s-t );
    
    XXis = XX( s, :, :); % 1 * p * Nk;
    XXis = permute(XXis, [3 2 1]);
    uuis = uis( s, :); % 1 * Nk
    uuis = repmat( uuis', [1 p]); % Nk * p
    Xuis = XXis .* uuis; % Nk * p
    Xuis = repmat( Xuis, [1 1 p]); % Nk * p * p
    
    XXit = XX( t, :, :); % 1 * p * Nk
    XXit = permute(XXit, [3 2 1]); % Nk * p
    uuit = uis( s, :); % 1 * Nk
    uuit = repmat( uuit', [1 p]); % Nk * p
    Xuit = XXit .* uuit; % Nk * p
    Xuit = permute( Xuit, [1 3 2]);
    Xuit = repmat( Xuit, [1 p 1]); % Nk * p * p
    
    XuuX_st1 = sum( Xuis .* Xuit, 1); % 1 * p * p
    XuuX_st(:,:, j) =  permute( XuuX_st1, [2 3 1]);
end
kern_mat = repmat( kern_st, [p p 1]);
kXuuX = sum( XuuX_st .* kern_mat, 3 )/(T * Nk);


vari = 1/(Nk * T) * pinv(Phi) * kXuuX * pinv(Phi);
end