function [num] = summary_RMSE_bias(B, B1_true)

Rep = size(B,2);
DEV = B - repmat( B1_true(:,1), [1, Rep]);
BIAS = mean( mean( abs( DEV ) ));
RMSE = mean( sqrt( mean( DEV.^2) ) );

num = [RMSE BIAS];
end
