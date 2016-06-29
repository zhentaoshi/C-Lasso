function [ post] = post_est_PLS_dynamic( T, g_data)
% Liangjun Su, Zhentao Shi, Peter Phillips (2016)
% collect the statistics after PLS for dynamic model

% INPUT:
%	T: number of time period
%	g_data: collect the data with group identity (group identity here is determined by CLasso)
%	class_a: coefficient estimated from CLasso

% OUTPUT: (each is a matrix containing the value of the coefficient, standard error and t-statistic.
%	classo: raw classo with bias correction
%	post: post lasso with bias correction

post_a = g_data.X \ g_data.y;

[bia] =  SPJ_PLS(T, g_data.y_raw, g_data.X_raw ) -post_a;
[vari] = var_PLS(T, post_a, g_data.y, g_data.X);
post.post_a_corr = post_a - bia;
post.se          = sqrt(diag(vari));
post.test_b      = post.post_a_corr./post.se ;
end
   