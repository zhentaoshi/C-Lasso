%% this is a script to generate the data for the classification


%% set the parameters
N = 100;
T = 40;

global p N_cut a0 K
p = 2;
N_cut =  N * [ 0.3, 0.6,  1 ];
K = length(N_cut);
a0 = [0.4, 1.6;      1, 1;    1.6, 0.4 ];

%% DGP1: the static linear regression model
[y, X] = DGP_static(N, T);
% y is an N*T vector
% X is an (N*T)*p matrix

% the first T observations belong to individual i = 1
% the following T observations belogn to individual i = 2
% and so on






% group0 = zeros(N, K);
% for k = 1:K
%     if k == 1
%         group0(1: N_cut(k) , k) = 1;
%     else
%         group0( (N_cut(k-1) + 1):N_cut(k), k ) = 1;
%     end
% end
% group0 = logical(group0);




% 
% post = zeros(K, p , Rep);
% post_var = post;
% cover   = post;



% 
% 
%     dat = dataset( kron( (1:N)', ones(T, 1) ), repmat((1:T)', N, 1), y, X );
%     dat.Properties.VarNames = {'N'  'T'  'y'  'X'};
% 
%     for k = 1:K
%         this_group = (group0(:, k)  == 1);
% 
%         index = 1:N;
%         g_index = index(this_group);
%         g_data = dat( ismember(dat.N, g_index), : );
% 		thisY = g_data.y;
% 		thisX = g_data.X;
% 
%         post_a = thisX\ thisY;
%         post(k, :, r) = post_a';
%         post_var(k,:,r) = sqrt( diag( var_PLS_static(T, post_a, y, X)));
%         cover(k,:,r) = abs( ( post_a' - a0(k,:) )./post_var(k,:,r) ) < 1.96;
%     end


% myTitle = [ 'PLS_static_Oracle_', num2str(N), '_T_', num2str(T), '_Rep_', num2str(Rep),'_seed_', num2str(seed), '.mat'];
% save( myTitle);


