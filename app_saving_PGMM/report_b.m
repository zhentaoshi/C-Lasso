function [beta_est, beta_est_co, group, group_co, unclass] = report_b( T, b, a, K )
% this script reports the classification outcomes.

% INPUT
%   b: ( N * p * K ) matrix
%   a0: (K * p) matrix
%   K: number of groups in the estimation

% OUTPUT
%   group is the uncoerced group.

global p tune_tol
 %% put each observation into a group. Coersion.

    N = size(b, 1);
    distance = zeros(N,K);
    
    tune_tol = 0.01 * p; % c_tol/( sqrt(N*T) * log( log(N*T) )); % the equality judgement follows the Liangjun's instruction

    group = zeros(N, K);
    group_co = zeros(N, K);
    for k = 1:K
        
        nominator = norms( bsxfun(@minus, b(:, :, k), a(k, :) ), 2 ,2 ); 
        % the nominator: inside the norm is (N * p). 
        % Outside the norm is ( N * 1 )
        distance(:,k) = nominator./norm( a(k,:) ) ;
         
    end
    [min_val, which] = min(distance, [], 2);
    
%%    % get group identity
    % judge if the minimum value is greater than the criterion.
    % if yes, classify; if no, do nothing and leave it blank.
    beta_est = 999 * ones(N, p );
    beta_est_co = 999 * ones(N, p );
    
    for i = 1:N
        group_co(i, which(i)) = 1;
        beta_est_co(i, :)  = a(which(i), :);
        if  min_val(i) < tune_tol
            group(i, which(i) ) = 1;
        end
    end
	group_co = logical(group_co);
    
%%
    index_n = 1:N;
    in_group = index_n( sum(group, 2) == 1 );
    out_of_group = setdiff( (1:N), in_group ); 
    unclass = length(out_of_group);
    

    %% assign to beta_est
    
    for i = in_group
        for k = 1:K
            if group(i,k) == 1
                beta_est(i,:) = a(k, :);
            end
        end
    end
    % for those unclassified K ~=1, average the unclassified observations.
    beta_est( out_of_group, : ) = mean( b( out_of_group, :, :), 3 );
end
