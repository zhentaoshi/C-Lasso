% this is a mimimum function to run PLS
load('PLS_data.mat') % load data

lambda_iter = 0.5* var(y) / nT^(1/3); % set the tuning parameter

% carry out PLS estimation
[b_iter, a_iter, group_iter] = SSP_PLS_est(nN, nT, y, X, K, lambda_iter, 80); 

