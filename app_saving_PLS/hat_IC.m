function h_Sigma  = hat_IC( y, X, b, a, K)
% Liangjun Su, Zhentao Shi, Peter Phillips (2014)
%
% this function computes "sigma_square_hat"
%
% INPUT
%   y: the dependent variable
%   X: the explantory variables
%   b: the estimate beta
%   a: the estimated alpha
%   K: number of groups
%
% OUPUT
%   h_Sigma: sigma_square_hat

N = size(b, 1);
T = size(y, 1)/N;

ds = dataset( kron( (1:N)', ones(T, 1) ), repmat((1:T)', N, 1), y, X );
ds.Properties.VarNames = {'N'  'T'  'y'  'X'};

% post-lasso estimation
if K == 1
    post_a = regress(ds.y, ds.X);
    post_b = repmat( post_a', N, 1) ;
else
    %% put each observation into a group.
    distance = zeros(N, K);
    for k = 1:K
        distance(:,k) = norms( bsxfun(@minus, b, a(k,:) ), 2, 2);
    end
    [ ~, group] = min(distance,[], 2); %
    %% post estimation
    post_a = zeros( size(a) );
    post_b = zeros( size(b) );
    
    for k = 1:K
        this_group = logical(group == k);
        
        index = 1:N;
        g_index = index(this_group);
        g_data = ds( ismember(ds.N, g_index), : ); % group-specific data
        
        
        post_a(k,:) = g_data.X \ g_data.y;
        
        
        post_b(this_group, :) = repmat(post_a(k, :), sum( this_group ), 1) ;
    end
end
%% IC with the post-selection estimator and penalization on the unclassified.
B = kron( post_b, ones(T,1)); % update the b as the post-estimator
h_Sigma =  mean( ( ds.y - sum( ds.X .*B, 2) ).^2 );
end