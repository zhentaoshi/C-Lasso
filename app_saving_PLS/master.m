% Liangjun Su, Zhentao Shi, Peter Phillips (2015)
% the master file of PLS estimation in the savings rate application

clear
global p K_max

cvx_solver mosek

IC_needed = 1;
tol = 0.0001;
R = 80;

load('balancedPanelX1995.mat')
X = [lagsaving, cpi, interest, gdp];
y = saving;


p = size(X, 2);
T = 15;
N = 56;

K_max = 5;
lamb.grid = 10;
lamb.min  = 0.2;
lamb.max  = 2.0;
lamb_const = lamb.min * (lamb.max / lamb.min ).^( ( (1:lamb.grid) - 1) /( lamb.grid -1 ) ); % the constant for lambda. very important!!
numlam = length(lamb_const);

index = dataset( code, year, y, X );
index.Properties.VarNames = {'N'  'T'  'y'  'X'};

y_raw = y;
X_raw = X;

for i = 1:N
    yi = y(index.N == i);
    mean_yi = mean(yi);
    yi = bsxfun(@minus, yi, mean(yi) );
    y(index.N == i) = yi/std(yi, 1);
    y_raw(index.N==i) = y(index.N == i) + mean_yi;
    
    Xi = X(index.N == i, : );
    mean_Xi = mean(Xi);
    Xi = bsxfun(@minus, Xi, mean(Xi) );
    X(index.N == i, :) = Xi./repmat( std(Xi, 1), [T 1] ) ;
    X_raw(index.N == i, :) = X(index.N == i, :) + repmat( mean(Xi), [T 1]);
end

ds = dataset( code, year, y, X, y_raw, X_raw );
ds.Properties.VarNames = {'N'  'T'  'y'  'X' 'y_raw' 'X_raw'};
%% initial values
beta_hat0 = zeros(N, p);
for i = 1:N
    yi = ds.y(ds.N == i );
    Xi = ds.X(ds.N == i, : );
    beta_hat0(i,:) = regress( yi , Xi );
end

%% estimation
TT = T;
IC_total = ones(K_max, numlam );

if IC_needed == 1
    for ll = 1:numlam
        disp(ll)
        
        a = ds.X \ ds.y; 
        bias = SPJ_PLS(T,ds.y_raw, ds.X_raw);
        a_corr = 2 * a - bias;
        IC_total(1, :) = mean( ( y - X*a_corr ).^2 );
        
        
        for K = 2:K_max
            Q = 999*zeros(K,1);
            
            lam = lamb_const(ll)*var(y) * T^(-1/3);
            [b_K, hat.a] = PLS_est(N, TT, y, X, beta_hat0, K, lam, R, tol); % estimation
            [~, H.b, ~, group] = report_b( b_K, hat.a, K );
            sum(group)            

            post_b = zeros(N, p);
            post_a = zeros(K, p);
            if K >=2
                for i = 1:K
                    NN = 1:N;
                    H.group = logical(group);
                    this_group = group(:,i);
                    if sum(this_group) > 0
                        g_index = NN(this_group);
                        g_data = ds( ismember(ds.N, g_index), : );

                        post = post_est_PLS_dynamic(T, g_data);
                        
                        e = g_data.y - g_data.X * post.post_a_corr ;
                        Q(i) = sum( e.^2 );
                        post_b(this_group,:) = repmat(post.post_a_corr', [sum(this_group), 1] );
                    end
                end
            end
            
            
            IC_total(K , ll) = sum(Q) / (N*T)
            
        end
    end
    %% calculate the IC
    pen = 2/3 * (N*T)^(-.5) * p .* repmat( (1:K_max)', [1 numlam]);
    IC_final = log(IC_total) + pen;
    disp(IC_final)
end

%% PLS estimation
K = 2;
lam = 1.5485 *var(y) * T^(-1/3);

[b_K, a] = PLS_est(N, T, y, X, beta_hat0, K, lam, R, tol);
[~, b, ~ , group] = report_b( b_K, a, K );

%% post estimation
est_lasso = zeros(p, 6);
est_post_lasso = zeros(p, 6);

for i = 1:K
    NN = 1:N;
    group = logical(group);
    this_group = group(:,i);
    g_index = NN(this_group);
    g_data = ds( ismember(ds.N, g_index), : ); % group-specific data
    post = post_est_PLS_dynamic(T, g_data);
    est_post_lasso(:,(3*i-2):(3*i)) =  [post.post_a_corr, post.se, post.test_b];
end

%% display the estimates
est_post_lasso = mat2dataset( est_post_lasso, 'VarNames', ...
    {'g1_coef', 'g1_sd', 'g1_t', 'g2_coef', 'g2_sd', 'g2_t'});
disp(est_post_lasso)

load('country56.mat')
country(group(:,1))
country(group(:,2))

g_PLS = zeros(56,1);
g_PLS( group(:,1) == 1 ) = 1;
g_PLS( group(:,2) == 1 ) = 2;

load('group_PGMM.mat')
g_PGMM = zeros(56,1);
g_PGMM( group_PGMM(:,2) == 1) = 2;
g_PGMM( group_PGMM(:,1) == 1) = 1;

sum(g_PLS == g_PGMM)
%% common FE

g_index = NN;
first_none_zero = min( NN );
g_data = ds( ismember(ds.N, g_index), : ); % group-specific data
post = post_est_PLS_dynamic(T, g_data);

[post.post_a_corr, post.se, post.test_b]

