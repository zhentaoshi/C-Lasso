% Liangjun Su, Zhentao Shi, Peter Phillips (2014)
%
% this script calculates the information criterion for the empirical
% application in the paper.
% the basic steps are similar to "master_est.m". Please refer to the
% comments in that file.

clear
global p

seed = 201;
rng(seed)

load('balancedPanelX1995.mat')
X = [lagsaving, cpi, interest, gdp];
y = saving;


p = size(X, 2);
T = 15;
N = 56;




% parameter for convergence
tol = 0.0001; % convergence tolerance level
R = 80; %  maximum number of iterations

%% demean and de-variance (not useful for simulation but useful for empirical application)
index = dataset( code, year, y, X );
index.Properties.VarNames = {'N'  'T'  'y'  'X'};

for i = 1:N
    yi = y(index.N == i);
    yi = bsxfun(@minus, yi, mean(yi) );
    y(index.N == i) = yi/std(yi, 1);
    
    Xi = X(index.N == i, : );
    Xi = bsxfun(@minus, Xi, mean(Xi) );
    X(index.N == i, :) = Xi./repmat( std(Xi, 1), [T 1] ) ; % update to the de-variance variables
end


% prepare the dataset. Useful for the functions.
ds = dataset( code, year, y, X );
ds.Properties.VarNames = {'N'  'T'  'y'  'X'};
%% initial values
beta_hat0 = zeros(N, p);
for i = 1:N
    yi = ds.y(ds.N == i );
    Xi = ds.X(ds.N == i, : );
    beta_hat0(i,:) = regress( yi , Xi ); %initial value
end


%%
TT = T;

K = 2;
lam = 1.5485*var(y) * T^(-1/3);

[b_K, a] = PLS_est(N, T, y, X, beta_hat0, K, lam, R, tol);
a



%% PLS estimation
Rep = 10;
optval = zeros(R, Rep);
QQ = zeros(1,Rep);
for rr = 1:Rep
    if rr == 1
    [~, a, optval_r] = PLS_est_Q(N, T, y, X, beta_hat0, beta_hat0, K, lam, R, tol);
    optval(:, rr) = optval_r;
    a
    else
    beta_random =  beta_hat0 + 2* rand(N,p)-1;
    [~, a, optval_r] = PLS_est_Q(N, T, y, X,  beta_random, beta_hat0 , K, lam, R, tol);
    end
    optval(:, rr) = optval_r;
    a
end

% export(dataset(group), 'file', 'group_id_PLS.csv', 'Delimiter', ',');
% b: b with coersion.
% group: group identity with coersion.
export(dataset(optval), 'file', 'opt_PLS.csv','Delimiter',',')

%% display the estimates






