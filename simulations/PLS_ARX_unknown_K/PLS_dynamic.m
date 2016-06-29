clear;
global p d D N_cut a0 K0

T = 80;
N = 100;
p = 1;





numlam = 2;
Rep = 2;


d = 2;
D = ( p - 1 ) + d;


N_cut =  N * [ 0.3, 0.6, 1 ];
K0 = length(N_cut);



SNR = 2;
a0 = SNR * [0.5 0.5 0.5; 0.4 0.4 0.4; 0.3, 0.3, 0.3 ];



tol = 0.0001;
R = 80;


lamb.grid = 5;
lamb.min  = 2;
lamb.max  = 5;
lambda = linspace( lamb.min, lamb.max, lamb.grid);

c.rho = [0.25:0.25:1, 1.5:0.5:3, 4:20];
rho = c.rho .* (log( N * T)/ (N * T) )  ;
K_max = 3;
K_grid = 3:K_max;



[g.K, g.rho, g.lambda] = meshgrid( K_grid, rho, lambda );
IC_data = dataset( g.lambda(:), g.K(:), g.rho(:));
IC_data.Properties.VarNames = {'lambda'  'K'  'rho'};
IC_data.IC1 = zeros(size(IC_data, 1), Rep);
IC_data.IC2 = zeros(size(IC_data, 1), Rep);


tic
for r = 1:Rep
    r
    toc
    [y, X, Dy, DX, Z] = DGP_dynamic(N, T);

    y = reshape(y, (T+d)*N, 1);

    X = reshape(X, (T+d)*N, p);

    for lam = lambda
        for K = K_grid
            [b_K, hat.a] = PLS_est(N, (T+d), y, X, K, lam, R, tol);
            hat.a
            [hat.b, hat.group, unclass] = report_b( b_K, hat.a, K, tol);

            [~, hat.IC, hat.post_a, hat.post_b] = hat_IC(y, X, hat.b, hat.a, unclass, K, rho);
            hat.post_a





            IC_data.IC1( (IC_data.lambda == lam) & (IC_data.K == K) , r ) = hat.IC(:,1);
            IC_data.IC2( (IC_data.lambda == lam) & (IC_data.K == K) , r ) = hat.IC(:,2);
        end
    end
end

toc
title = ['IC_AR_x_',num2str(N),'_T_', num2str(T),'_rep_', num2str(Rep)];
save(  [title, '.mat'], 'IC_data');
export( IC_data, 'file',  [title, '.csv'],'Delimiter', ',') ;


