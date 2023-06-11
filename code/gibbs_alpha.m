%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Article : Multiclass classification of growth curves using random change 
%            points and heterogeneous random effects
%
%  Authors : Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and
%            Scott A. Sisson
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gibbs sampling for mean of normal distribution describing random
% intercept
function output = gibbs_alpha(ID, y, xf, xr, beta_i, r_int, precision, ...
                              same_indv, prior_alpha)

A0       = prior_alpha.mat;
a0       = prior_alpha.mean;
num_pred = size(xf, 2);
x_var_x  = zeros(num_pred, num_pred);
x_var_y  = zeros(num_pred, 1);

for i=1:ID(end)
    index   = same_indv{i};
    x_var_x = x_var_x + xf(index,:)' * (precision(i)) * xf(index,:);
    x_var_y = x_var_y + xf(index,:)' * (precision(i)) * (y(index) - ...
              xr(index,:) * beta_i(:,i) - r_int(i));
end

inv_A0       = inv(A0);
post_cov_mat = inv(x_var_x + inv_A0);
post_mean    = post_cov_mat * (x_var_y + inv_A0 * a0);
output       = mvnrnd(post_mean, post_cov_mat);
end