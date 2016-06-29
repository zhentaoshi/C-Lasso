% Su, Shi, Phillips (2015)
% this script applies PLS on the data provided by Bonhomme and Manresa
% (2015) about democracy.

clear;
load('data84.mat')

X = [ldem linc];
y = dem;

IC_needed = 1;

[mean(y), median(y), std(y), min(y), max(y);
    mean(linc), median(linc), std(linc), min(linc), max(linc)];

xx = reshape(X(:,1), [7 84]);
country_name = unique(country);

T = 7;
N = length(y)/T;

code = kron( (1:N)', ones(T,1));
year = repmat( (1:T)', [N 1]);
%%
global p R tol
p = size(X, 2);

lamb.grid = 10;
lamb.min  = 0.2;
lamb.max  = 2;
lamb_const = lamb.min * (lamb.max / lamb.min ).^( ( (1:lamb.grid) - 1) /( lamb.grid -1 ) );
lambda = lamb_const* T^(-1/3);
numlam = length(lambda);
disp(lamb_const);

K_max = 5;
tol = 0.0001;
R = 80;

%% 

y_raw = y;
X_raw = X;

index = dataset( code, year, y, X );
index.Properties.VarNames = {'N'  'T'  'y'  'X'};

for i = 1:N
    yi = y(index.N == i);
    yi = demean(yi);
    y(index.N ==i ) = yi;
    
    Xi = X(index.N == i, : );
    Xi = demean(Xi);
    X(index.N == i, :) = Xi;
end

y = y/std(y);
y_raw = y_raw/std(y);
X(:,1) = X(:,1)/std(X(:,1));
X_raw(:,1) = X_raw(:,1) / std(X(:,1));
X(:,2) = X(:,2)/std(X(:,2));
X_raw(:,2) = X_raw(:,2) / std(X(:,2));




ds = dataset( code, year, y, X, y_raw, X_raw );
ds.Properties.VarNames = {'N'  'T'  'y'  'X' 'y_raw' 'X_raw'};

%% initial values
beta_hat0 = repmat( regress(y,X)', [N 1]);
for i = 1:N
    yi = ds.y(ds.N == i );
    Xi = ds.X(ds.N == i, : );
    if ( var(yi) > 0.001)
        beta_hat0(i,:) = regress( yi , Xi );
    end
end

TT = T;
IC_total = ones(K_max, length(lambda) );
lambda = lambda * var(y);

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
%% estimation

lam = 1.1990;
K = 3;

[b_K, a] = PLS_est(N, T, y, X, beta_hat0, K, lam, R, tol);
[~, b, ~ , group] = report_b(b_K, a, K );

est_post_lasso = zeros(p, 6);

for i = 1:K
    NN = 1:N;
    group = logical(group);
    this_group = group(:,i);
    g_index = NN(this_group);
    g_data = ds( ismember(ds.N, g_index), : );
    [ post] = post_est_PLS_dynamic(T, g_data);
    
    est_post_lasso(:,(3*i-2):(3*i)) =  [post.post_a_corr, post.se, post.test_b];
end

for kk = 1:K
    disp(country_name( group(:, kk) ) )
end
sum(group)


%% common FE
g_index = NN;
g_data = ds( ismember(ds.N, g_index), : ); 
[post] = post_est_PLS_dynamic(T, g_data);

common =  [post.post_a_corr, post.se, post.test_b];

out =  [common, est_post_lasso]

