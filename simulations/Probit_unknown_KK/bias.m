function [bia] = bias(a, y, X, Z, MT)










global T
n = size(y, 1)/T


ZZ =
XX =

Q = XX * ZZ';
A = 1/(N * T) * Q' * Q;


u = y - X * a;

XZZU = @(s, t) X(s) * Z (s) * Z(t) * u(t);
kernel = @( u ) (1- abs(u) / MT ) * ( abs(u) <= MT );
comb_st = combination(1:T, 1:T);


num2cell(comb_st)
kXZZU = cellfunc( @(s,t) kernel( s - t ) * XZZU(s, t), cell );
B = 1/(n^0.5 * T^(1.5) ) * mean(kXZZU);

bia = 1/( n * T)^0.5 * pinv(A) * B;

end
