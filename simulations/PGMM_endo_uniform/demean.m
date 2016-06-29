function [y] = demean(y)


    T = size(y, 1);
    if ( size(y, 2) == 1 )
        y = y - mean(y);
    else
        y = y - repmat( mean(y), T, 1);
    end
end
