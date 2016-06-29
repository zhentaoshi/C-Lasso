function y = delta(x)



global T

x1 = x( 1:(T-1), : );
x2 = x( 2:T    , : );

y = x2 - x1;
end
