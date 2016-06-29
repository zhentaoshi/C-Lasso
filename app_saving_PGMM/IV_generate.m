function [Dy, DX, Z] = IV_generate(y, X )
% this function generates a T-period AR(1) sequence for each individual i.
% It accepts exogenous regressors

% because of the dynamics, 
% consume one period in differentiation
% consume two more periods for over-identifying IV
    
    
    % OUPUT
    %   y: the dependent variable. In our context it is the first diff. of
    %   the orignal time series.
    %   X: the explanatory variables. In our context it is the lag term of
    %   the dependent variable as well as possible exogenous regressors.
    %   Z: instruments.
    
    global  p N T d
    
    % because of the lag, it will indeed waste 1 + d observations. 
    % we lose 2 because of the Dy and Dy_lag
    % we lose (d-1) from the IVs.
    
    y_mat = reshape(y, [T N]);
    X_mat = reshape(X, [T N p]);
    
    D_y_mat = y_mat( 2:T, :) - y_mat(1:(T-1), :);
    % but the diff in X still eat up one period
    D_X_mat = X_mat( 2:T, :, :) - X_mat(1:(T-1), :, :);
    
    % use level as instrument
    IV_y1 = y_mat(1:(T-2), :);
    IV_y1 = [zeros(1,N); IV_y1];
    IV_y2 = X_mat(1:(T-2), :, 1);
    IV_y2 = [zeros(1,N);IV_y2];
    
    Delta_y = IV_y1 - IV_y2;
    std_Delta_y = std( reshape(Delta_y, [], 1) );
    Z(:,:,1) = IV_y1 -IV_y2;
    Z(:,:,2) = std_Delta_y; % to keep the two IVs of the same order

    Z(:,:,d+1:(p+d-1)) = D_X_mat(:, :, 2:p);


    % clean the zeros
    D_y_mat(1, :) = [];
    D_X_mat(1, :,:) = [];
    Z(1, :, :) = [];
   
    Dy = D_y_mat;
    DX = D_X_mat;
    
   % Z(:, :, 2) = [];
end
    
    
    
    


