clear;
global N T p N_cut a0 sigma_x sigma_e K_max

T = 20;
N = 200;
p = 2;

N_cut =  N * [ 0.3, 0.6,  1 ];



SNR = 2;
a0 = SNR * [0.25, 0.75; 0.5, 0.5; 0.75, 0.25];

sigma_x = 1; sigma_e = 1;


tol = 0.0001;
R = 80;


lamb.grid = 10;
lamb.min  = 0.1;
lamb.max  = 20;
lambda = lamb.min * (lamb.max / lamb.min ).^( ( (1:lamb.grid) - 1) /( lamb.grid -1 ) );
c.rho = 1:50;
rho = c.rho .* (log( N * T)/ (N * T) )  ;
K_max = 6;
K_grid = 1:K_max;


Rep = 2;
[d.K, d.rho, d.lambda] = meshgrid( K_grid, rho, lambda );
IC_data = dataset( d.lambda(:), d.K(:), d.rho(:));
IC_data.Properties.VarNames = {'lambda'  'K'  'rho'};
IC_data.IC = zeros(size(IC_data, 1), Rep);


for r = 1:Rep
    r
[y, X] = DGP_static();

    for lam = lambda
        for K = K_grid

            [b_K, hat.a] = PLS_est(y, X, K, lam, R, tol);
            [hat.b, hat.group, unclass] = report_b( b_K, hat.a, K);

            [hat.Sigma, hat.IC] = hat_IC(y, X, hat.b, hat.a, unclass, K, rho);


            IC_data.IC( (IC_data.lambda == lam) & (IC_data.K == K) , r ) = hat.IC;
        end
    end
end

save( ['IC_static_N_',num2str(N),'_T_', num2str(T),'_rep_', num2str(Rep), '.mat'] , 'IC_data');
export( IC_data, 'file',  ['IC_static_no_penalty_unclass_N_',num2str(N),'_T_', num2str(T),'_rep_', num2str(Rep), '.csv'],...
   'Delimiter', ',') ;



