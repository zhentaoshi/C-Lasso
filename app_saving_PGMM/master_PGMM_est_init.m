% the paths produced by various initial values

clear all
global p T N d ;
cvx_solver mosek 


% load the data
load('balancedPanelX1995.mat')
X = [lagsaving, cpi, interest, gdp];
y = saving;

p = size(X, 2);
T = 15;
N = 56;
d = 2; % number of IVs

seed = 501;
rng(seed);


% parameter for convergence
tol = 0.001; % convergence tolerance level
R = 80; %  maximum number of iterations
W = 1; % GMM weighting matrix. Equal weighting.

%% scale the dependent variable and its lag term
[Dy, DX, Z ] = IV_generate(y, X);
std_yi = std(Dy);
Dy = bsxfun(@rdivide, Dy, std_yi) ;
DX = bsxfun(@rdivide, DX, std(DX, 1));
Z(:,:,3:5) = DX(:,:,2:4);

K = 2;

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


lam = 0.7188 * T^(-1/3)* var(reshape(Dy, [],1 ));




%% initial value

Rep = 10;
optval = zeros(R, Rep);
QQ = zeros(1,Rep);
for r = 1:Rep
    if r == 1
       [~, hat.a, optval(:,r) ] = PGMM_est_Q(Dy, DX, Z,  beta_hat0, W, K, lam, R, tol);
       hat.a
    else
        beta_random =  beta_hat0 +   2 * rand(N,p) - 1;
        [~, hat.a, optval(:,r)] = PGMM_est_Q(Dy, DX, Z,  beta_random, W, K, lam, R, tol);
    end
    hat.a
end

export(dataset(optval), 'file', 'opt_PGMM.csv','Delimiter',',')

