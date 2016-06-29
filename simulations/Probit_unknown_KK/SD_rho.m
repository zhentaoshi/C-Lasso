function sd_rho = SD_rho( hat_a, fe, N, T, y, X)









YL = y; YR = X(:,:,1); X = X(:,:,2);

[logl, grad, Hess,  Grho, Gbeta]=G_H(hat_a, fe,YL,YR,X);


Grho_vector = reshape(Grho, [T*N 1]);
Gbeta_vector= reshape(Gbeta,[T*N 1]);

Omega = cov([Grho_vector, Gbeta_vector]);

vari =   pinv(Hess) * Omega * pinv(Hess) ;

sd_rho = sqrt( vari(1,1)/(N*T) );


end



