global p N_cut a0
load('IC_PNL_dyn_200_T_25_rep_500_seed_22000')
rho1 = 1/4 * log(log(T))* T^(-1) * p;
pen = rho1 * repmat(1:5, [r 1]);
IC_pen = IC_data + pen;
[ ~, II] = min(IC_pen,[], 2);

STAT = zeros(5,1);
for kk = 1:5
    STAT(kk) = mean(II == kk);
end

load('PNL_real_dyn_200_T_25_rep_500_seed_22000.mat')



workhorse
