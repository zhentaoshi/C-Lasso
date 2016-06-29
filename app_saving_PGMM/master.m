% Su, Shi and Phillips (2015)
% PGMM estimation of the savings rate application

clear
global p T N d pi;
IC_needed = 1;
cvx_solver mosek
load('balancedPanelX1995.mat')
X = [lagsaving, cpi, interest, gdp];
y = saving;

p = size(X, 2);
T = 15;
N = 56;
d = 2; % number of IVs

K_max = 5;
lamb.grid = 10;
lamb.min  = .2;
lamb.max  = 2.0;
lamb_const = lamb.min * (lamb.max / lamb.min ).^( ( (1:lamb.grid) - 1) /( lamb.grid -1 ) ); % the constant for lambda. very important!!

lambda = lamb_const* T^(-1/3);
numlam = length(lambda);

tol = 0.0001;
R = 80;
W = 1;

%% scale the dependent variable and its lag term
[Dy, DX, Z ] = IV_generate(y, X);
std_yi = std(Dy);
Dy = bsxfun(@rdivide, Dy, std_yi) ;
DX = bsxfun(@rdivide, DX, std(DX, 1));
Z(:,:,3:5) = DX(:,:,2:4);

beta_hat0 = zeros(N, p);
for i = 1:N
    DXi = permute( DX(:,i,:), [1 3 2]);
    Zi = permute( Z(:,i,:), [1 3 2]);
    Dyi = Dy(:, i);
    beta_hat0(i,:) = pinv(DXi' * (Zi * Zi') * DXi) * (DXi' * (Zi * Zi') * Dyi); %initial value
end
beta_hat00 = beta_hat0;

T = 13;
ds = dataset( kron( (1:N)', ones(T,1) ), repmat( (1:T)', [N 1]),...
    reshape(Dy, [N*T 1]), reshape(DX, [N*T size(DX,3)] ), reshape(Z, [N*T size(Z,3)]) );
ds.Properties.VarNames = {'N'  'T'  'y'  'X', 'Z'};


%%
if IC_needed == 1
    IC_value = zeros(K_max, numlam);
    IC = IC_PGMM(ds, ones(N,1)); % K = 1 case
    IC_value(1, :) = log( IC/(N*T) );
    
    for ll = 1:numlam
        disp(ll)
        for K =  2:K_max
            lam = lambda(ll)*var(ds.y);
            [b_K, hat.a] = PGMM_est(Dy, DX, Z, beta_hat0, W, K, lam, R, tol);
            [~, hat.b, ~ , group] = report_b( T, b_K, hat.a, K );
            sum(group)
            IC_kk = zeros(1,K);
            for  kk = 1:K
                IC_kk(kk) = IC_PGMM(ds, group(:,kk) );
            end
            IC_value(K,ll) = log( sum(IC_kk)/(N*T) ) 
        end
    end
    
    
    pen =  (2/3)*(N*T)^(-.5)* p*(1:K_max);
    IC_final = IC_value + repmat( pen, [numlam, 1])';
    [~, l_hat] = min(IC_final(2,:));
end

%% the information criterion selects
% K = 2 and lambda(8)

K = 2;
l_hat
lam = lambda(l_hat)*var(ds.y);


[b_K, hat.a] = PGMM_est(Dy, DX, Z,  beta_hat0, W, K, lam, R, tol);
[~, hat.b, ~ , group] = report_b( T, b_K, hat.a, K );
H = hat_IC_PGMM( ds, hat.b, hat.a, K, group);
sum(group)


est_post_lasso = zeros(p,6);
pi =  0.05;
for kk = 1:K
    this_group = group(:,kk);
    dat = ds;
    a_hat = hat.a(kk,:)';
    
    index = 1:N;
    g_index = index(this_group);
    g_data = dat( ismember(dat.N, g_index), : ); % group-specific data
    
    ky = g_data.y;
    kX = g_data.X;
    kZ = g_data.Z;

    %%
    Nk = sum(this_group);
    XZZX = (kX' * kZ) * pinv(kZ' * kZ) * (kZ' * kX)/(Nk*T)+ridge(pi, N, T);
    XZZy = (kX' * kZ) * pinv(kZ' * kZ) * (kZ' * ky)/(Nk*T);

    post_a = (pinv(XZZX) * XZZy)';
    
    [vari] = var_post_GMM(T, post_a', ky, kX, kZ, pi); % no ridge used
    se = sqrt(diag(vari));
    est_post_lasso( :, (3 * ( kk-1 ) + 1) :3 * kk ) = [post_a; se'; post_a./se' ]';
    
end


est_post_lasso = mat2dataset( est_post_lasso, 'VarNames', ...
    {'g1_coef', 'g1_sd', 'g1_t', 'g2_coef', 'g2_sd', 'g2_t'});
disp(est_post_lasso)

load('country56.mat')
country(group(:,1))
country(group(:,2))
sum(group)

%% pooled estimation
% when the two groups are put together, we found the IV's is weak
% this is the heterogeneity induced weak IV.
% we add a positive value to the matrix to stablize the inverse, as in
% ridge regression

XZZX = ds.X' * ds.Z * pinv(ds.Z' * ds.Z) * ds.Z' * ds.X/(N*T) + ridge(pi, N, T);
XZZy = ds.X' * ds.Z * pinv(ds.Z' * ds.Z) * ds.Z' * ds.y/(N*T) ;

alpha1 = XZZX \ XZZy;

[vari] = var_post_GMM(T, alpha1, ds.y, ds.X, ds.Z,pi);
se = sqrt(diag(vari));
[alpha1, se]

group_PGMM = group;
save('group_PGMM.mat', 'group_PGMM');
