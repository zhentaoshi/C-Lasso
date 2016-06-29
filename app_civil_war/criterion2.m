function [d] = criterion2( a_old, a_new, optval_old, optval_new, tol)
% Liangjun Su, Zhentao Shi, Peter Phillips (2014)
% the convergence criterion of the algorithm

% INPUT
%   a_old: a in the previous round
%   a_new: a in the current round
%   tol: tolerenace level

% OUPUT
%   d: a dummy variable. d = 1 indicates "converged". Otherise d = 0

    d = 0;
    
    a_nominator = sum( norms( abs( a_old - a_new ), 2,2 ) );
    a_denominator = sum( norms( abs( a_old ), 2, 2 ) ) + 0.001;

    %%
    if a_nominator/ a_denominator < tol &&  (optval_old - optval_new < tol)
        d = 1;
    end

end