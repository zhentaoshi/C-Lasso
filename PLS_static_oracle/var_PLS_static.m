function [vari] = var_PLS_static(T, a, y, X )


Nk = size(y, 1)/T;

Phi = 1/(Nk * T) * (X' * X);

ui = y - X*a;
ui = demean(ui);
Xu = X .* repmat(ui, [1 size(X,2)]);
Omega = Xu' * Xu/(Nk * (T-1) );

vari = 1/(Nk * T ) * pinv(Phi) * Omega * pinv(Phi);
end
