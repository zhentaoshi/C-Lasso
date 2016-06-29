function [IC,alpha_hat] = IC_PGMM(ds, this_group)

global  T
if sum(this_group) == 0;
    IC = 0;
else
    sel = logical( kron( this_group, ones(T,1) ));
    n = sum(sel)/T;
    p = size(ds.X,2);

    y = ds.y(sel, :);
    X = ds.X(sel, :);
    Z = ds.Z(sel, :);
    
    cvx_begin
        variable a(p,1);
        g = Z' * (y - X * a);
        Q = quad_form(g, 1 );
        
        minimize(Q)
        subject to 
         %  a(1) <= .95;
         %  a(1) >= 0;
            
        cvx_end

    alpha_hat = a;

    
    e = y - X * alpha_hat;
    % eZ = bsxfun(@times, e, Z);
    % IC = trace(eZ' * eZ);
    IC = e'*e;
end