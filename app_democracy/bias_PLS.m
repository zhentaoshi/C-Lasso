function [bia] = bias_PLS(T, a, y, X )
% Liangjun Su, Zhentao Shi, Peter Phillips (2014)
% 
%%
% this function calculate A and B, the two components of the bias term as
% in the note p.2.
%
% INPUT
%   y, X, Z are not those from the full dataset. They are from one group.
%   a: the post-estimator without bias correction.
%
% OUTPUT
%   bia: bias to be corrected
%   vari: the variance of the bias-corrected estimator

global p
MT = ceil( 2* T^(1/4)); % use 1/5 or 1/8.
Nk = size(y, 1)/T; % the number of witin-group observaations

%% 
Phi = 1/(Nk * T) * (X' * X); % correct

Xa = X * a;
ui = y - Xa; 
%% 
uis = zeros( T, Nk);
XX = zeros( T, p, Nk);

for i = 1:Nk
    u_within_group = ui( ( T*(i-1) + 1 ): (T*i) );
    uis( :, i) = bsxfun(@minus, u_within_group, mean(u_within_group));
    XX(:, :, i) = X( ( T*(i-1) + 1 ): (T*i), : );   
end

%%
kernel = @( u ) (1- abs(u) / MT ) .*  ( abs(u) <= MT );
combNK = combnk(1:T, 2);
comb_st = [ repmat( (1:T)', 1, 2  ); combNK; combNK( :, [2,1] )  ]; % all combinations, including self

TTT = size(comb_st, 1);
kern_st = zeros(1, TTT);
XU_st = zeros(p, TTT);

for j = 1:TTT
    s = comb_st(j,1);
    t = comb_st(j,2);
    kern_st(j) = kernel( s-t );
    
    XXis = XX( s, :, :); % 1 * p * Nk
    XXis = permute(XXis, [ 3 2 1]); % Nk * p
    uuit = uis( t, :); % 1 * Nk
    
    XU_st(:, j) =   XXis' * uuit';
end

kXU = sum( bsxfun(@times, XU_st, kern_st), 2);
B = -1/(Nk^0.5 * T^(1.5) ) * kXU; 
bia = (1/(Nk * T)^0.5) * pinv(Phi) * B;
end