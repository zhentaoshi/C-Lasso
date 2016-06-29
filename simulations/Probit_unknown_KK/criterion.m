function [d] = criterion( a_old, a_new, b_old, b_new, tol)












    d = 0;
    a_nominator = sum( abs( a_old - a_new ) );
    a_denominator = sum( abs( a_old ) ) + 0.001;

    b_nominator = mean( abs( b_old - b_new ) );
    b_denominator = mean( abs( b_old ) ) + 0.001;

    b_nominator / b_denominator;


    if a_nominator/ a_denominator < tol &&  b_nominator / b_denominator < tol
        d = 1;
    end

end
