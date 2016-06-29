
global p N_cut a0 sigma_x sigma_e K_max static

tic
static = 1;

N = 100;
T = 10;
Rep = 2;
numlam = 5;
seed = 1;
rng(seed);


p = 2;

N_cut =  N * [ 0.3, 0.6,  1 ];



SNR = 2;
a0 = SNR * [0.25, 0.75; 0.5, 0.5; 0.75, 0.25];

sigma_x = 1; sigma_e = 1;


tol = 0.0001;
R = 80;


lamb.grid = numlam;
lamb.min  = 0.1;
lamb.max  = 10;

lambda = lamb.min * (lamb.max / lamb.min ).^( ( (1:lamb.grid) - 1) /( lamb.grid -1 ) );


c.rho = [0.25:0.25:1, 1.5:0.5:3, 4:20];
rho = c.rho .* (log( N * T)/ (N * T) )  ;
K_max = 3;
K_grid = 1:K_max;



[d.K, d.rho, d.lambda] = meshgrid( K_grid, rho, lambda );
IC_data = dataset( d.lambda(:), d.K(:), d.rho(:));
IC_data.Properties.VarNames = {'lambda'  'K'  'rho'};
IC_data.IC1 = zeros(size(IC_data, 1), Rep);
IC_data.IC2 = zeros(size(IC_data, 1), Rep);



 for r = 1:Rep

    [y, X] = DGP_NL(N, T);



    beta_hat0 = zeros(N,p);






     y = reshape(y, T*N, 1);
     X = reshape(X, T*N, p);

    for lam = lambda
        disp(lam)
        for K = K_grid

            [b_K, a_out, c_out] = PNL_est(N, T, beta_hat0, y,  X, K, lam, R, tol );
            hat.a = a_out;
            hat.c = c_out;
            disp(a_out);
            [hat.b, hat.group, unclass] = report_b( b_K, hat.a, K, tol);

            [hat.Sigma, hat.IC, group,~, ~] = hat_IC_NL(N, y, X, hat.b, hat.a, hat.c, unclass, K, rho);





            IC_data.IC1( (IC_data.lambda == lam) & (IC_data.K == K) , r ) = hat.IC(:,1);
        end
     end

	toc
	if mod(r,min(100, Rep) ) == 0
		title = ['IC_PNL_',num2str(N),'_T_', num2str(T),'_rep_', num2str(r),...
			'_seed_', num2str(seed) ];
		save(  [title, '.mat'], 'IC_data');

	end
 end

