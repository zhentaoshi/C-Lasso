clear
global p a0 K_max

tol = 0.0001;
R = 50;

%%
load('input_balanced_102.mat')
ds1 = dataset('XLSFile','country_name_code.xlsx');
p = 3;

%% remove the all-one all-zero countries
country102 = ds1.country( ismember(ds1.code, codelist) );
country102 = country102(  mean(y) < .99 &  mean(y) >= 0.01 );
code102 = ds1.code( ismember(ds1.code, codelist) );
code102 = code102( mean(y) <.99 & mean(y) >=0.01 );
wanted = ismember(codelist, code102);

T = 39; 
% 1999 has no information of population
% we use data 1960-1999
rho = 1/3 * log(log(T))* T^(-1) * p;

y = y(:, wanted);
var_y = var(y);
ylag = y(1:(T-1), :);
y    = y(2:T, :);

N = sum(wanted);

%% summary statistic
lpop =  ( x1(2:(T), wanted) )- (x1(1:(T-1), wanted));
gdp = log( x2(2:(T), wanted))- log( x2(1:(T-1), wanted));

lpop_vector = reshape( lpop, [N*(T-1) 1] );
gdp_vector = reshape(  gdp, [N*(T-1) 1]);
y_vector = reshape(  y, [N*(T-1) 1]);
ylag_vector = reshape(  ylag, [N*(T-1) 1]);

d0 = [y_vector ylag_vector  gdp_vector lpop_vector];
summary =[ mean(d0); median(d0); std(d0); min(d0); max(d0) ]';

%% scale normalization
sd_y = std(y_vector);

lpop  = lpop./repmat(std(lpop), [T-1, 1])*sd_y;
gdp  = gdp./repmat(std(gdp), [T-1, 1])*sd_y;

X = zeros(T-1, N, p);
X(:,:,1) = ylag;
X(:,:,2) = demean(gdp);
X(:,:,3) = demean(lpop);

%% tuning parameter
K_max = 2;
lamb.grid = 10;
lamb.min  = 0.01;
lamb.max  = 0.1;
lamb.const = lamb.min * (lamb.max / lamb.min ).^( ( (1:lamb.grid) - 1) /( lamb.grid -1 ) ); % the constant for lambda. very important!!
numlam = length(lamb.const);

T = T - 1;
K_grid = 2:K_max;
a0 = zeros(K_max, p); % is not useful, but to make a script getting through. No effect.
HAT_IC = zeros(K_max, numlam);
%%

X_all= reshape(X, T*N, p);
y_all = reshape(y, T*N,1);
beta_temp = fitglm( X_all , y_all, 'distr', 'binomial', 'link', 'probit' );
beta_temp_coef = beta_temp.Coefficients.Estimate;
y_star_pred = [ X_all, ones(size(X_all,1),1)] * beta_temp_coef;
y_pred = (y_star_pred >= 0);
pseudo_R = mean(y_pred .* y_all);
sign_y_all = (2*y_all - 1);


%% initial value
%     in each group we run a probit regression

% we tried to start from the QMLE initial, but
% some values are extraordinary, for example, those countries with 
% very few civil war incidence.
% we instead start from the pooled FE probit.

%% K = 1 case
y = reshape(y, T*N, 1);
X = reshape(X, T*N, p);

cvx_begin
variable c(N, 1);
variable b(p,1)

B1 = kron( b, ones(T,1) );
cst = kron(c, ones(T,1) );
Q =  -1/(N*T) * sum( log_normcdf( ( cst + X * b ) .* sign_y_all ) );
minimize( Q );

cvx_end
HAT_IC(1,:) = 2*Q;

%% CLasso
alpha_hat0 = zeros(N,1);
beta_hat0 = repmat(b', [N 1]);

for ll = 1:numlam
    lam = lamb.const(ll) * var(y) * T^(-1/3);
    for K = K_grid
        [b_K, a_out, c_out] = PNL_est(N, T, beta_hat0, y,  X, K, lam, R, tol ); % estimation
        hat.a = a_out;
        hat.c = c_out;
        [hat.b, hat.group] = report( b_K, hat.a, K);
        
        loglik = hat_IC_NL(N, y, X, hat.group, hat.a, hat.b, K);
        HAT_IC(K, ll) = 2 * loglik;
    end
end

%% report information criterion
pen = rho * repmat( ((1:K_max)-1)', [1 numlam]);
HAT_IC1 = HAT_IC + pen;
disp(HAT_IC1)
%% the information criterion pick K = 2 and lambda = lambda.const(7)

lam = lamb.const(7) * var(y) * T^(-1/3);
K = 2;

[b_K, a_out, c_out] = PNL_est(N, T, beta_hat0, y,  X, K, lam, R, tol ); % estimation
hat.a = a_out;
hat.c = c_out;
[hat.b, hat.group] = report( b_K, hat.a, K);
[loglik, ~, group, ~, post_b] = hat_IC_NL(N, y, X, hat.group, hat.a, hat.b, K);

%% report membership
y_mat = reshape(y, [T N]);
mean_y_mat = mean(y_mat);
for k = 1:K
    dataset( [mean_y_mat(hat.group == k)]', country102( hat.group == k) )
    code_k_long = logical( kron( (hat.group == k), ones(T,1)) );
    y_this_group = y(code_k_long);
    [ mean( y_this_group ), std( y_this_group )/sqrt(length(y_this_group)) ]
end

%% report coefficients

b % FE probit model one group
hat.a % FE probit model with C-Lasso

%%
X_mat(:,:,1) = reshape( X(:,1), [T N]);
X_mat(:,:,2) = reshape( X(:,2), [T N]);
X_mat(:,:,3) = reshape( X(:,3), [T N]);

coef = zeros(3,6);
disp(beta_temp)

% common estimator

this_group = logical(ones(size(hat.group)));
this.y = y_mat(:, this_group);
this.X = X_mat(:, this_group, :);
Nk = sum(this_group);

y_vector = reshape(this.y, T*Nk, 1);
X_vector = reshape(this.X, T*Nk, p);
sign_this_y = 2*y_vector - 1;


[post_a, post_c] = solve( zeros(3,1), Nk, T, X_vector, y_vector);
SD_a1 = SD_rho_robust(post_a, post_c, Nk, T, this.y, this.X);
bias1 =  SPJ( post_a, y_vector, X_vector, Nk, T)- post_a;

POST_a_corr = post_a - bias1;
coef(:,1:2) = [POST_a_corr, SD_a1];
[POST_a_corr, SD_a1, POST_a_corr./ SD_a1]


for kk = 1:K
    this_group = logical (hat.group == kk);
    this.y = y_mat(:, this_group);
    this.X = X_mat(:, this_group, :);
    Nk = sum(this_group);
    
    y_vector = reshape(this.y, T*Nk, 1);
    X_vector = reshape(this.X, T*Nk, p);
    X_vector(:, 2:3) = X_vector(:,2:3);
    sign_this_y = 2*y_vector - 1;    
    
    [post_a, post_c] = solve( zeros(3,1), Nk, T, X_vector, y_vector);
    
    SD_a1 = SD_rho_robust(post_a, post_c, Nk, T, this.y, this.X);
    bias1 =  SPJ( post_a, y_vector, X_vector, Nk, T)- post_a;
    
    POST_a_corr = post_a - bias1;
    coef(:, (2*kk+1):(2*kk+2)) = [POST_a_corr, SD_a1];
    [POST_a_corr, SD_a1, POST_a_corr./SD_a1]
end;


