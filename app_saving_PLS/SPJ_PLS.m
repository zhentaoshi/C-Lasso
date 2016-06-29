function [out, vari] = SPJ_PLS(T, y_vector, X_vector)
% half panel Jackknife
% develop on Dhaene and Jochmans (2015)

p = size(X_vector,2);
Nk = size(y_vector,1) / T;

period1_i = [ones(floor(T/2), 1); zeros( ceil(T/2),1)];
period1 = logical( kron(  ones(Nk, 1), period1_i )) ;
period2 = logical( 1-period1);

theta_bar = zeros(p, 2);
V = zeros(p,2);

for tt = 1:2
    if tt == 1
        X_half = X_vector(period1, :);
        y_half = y_vector(period1);
    elseif tt == 2;
        X_half = X_vector(period2, :);
        y_half = y_vector(period2);
    end
    half_t = size(y_half, 1)/Nk;
    y_mat = reshape(y_half, [half_t Nk]);
    X_mat = reshape(X_half, [half_t Nk p]);
    
    y_mat = demean(y_mat);
    for pp = 1:p
        X_mat(:,:,pp) = demean( X_mat(:,:,pp) );
    end
    
    y_half = reshape( y_mat, [half_t*Nk  1]);
    X_half = reshape( X_mat, [half_t*Nk  p]);
    
    b = X_half\y_half;
    theta_bar(:,tt) = b;
end
out =  mean(theta_bar, 2);
end