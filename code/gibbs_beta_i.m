%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Article : Multiclass classification of growth curves using random change 
%            points and heterogeneous random effects
%
%  Authors : Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and
%            Scott A. Sisson
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gibbs sampling for random slopes
function output = gibbs_beta_i(ID, y, xf_alpha, xr, r_int, precision, ...
                               grp, grp_mean, grp_cov_mat, same_indv)             

num_pred         = size(xr, 2);
y_minus_xf_alpha = y - xf_alpha;
output           = zeros(num_pred, ID(end));
zero_mat_x       = zeros(num_pred, num_pred);
zero_mat_y       = zeros(num_pred, 1);
i_mat            = eye(num_pred);
prior_precision  = grp_cov_mat;
prior_mean       = grp_mean;
active_grp       = unique(grp);
random           = mvnrnd(zero_mat_y', i_mat, ID(end));

for i=1:length(active_grp)
    j                      = active_grp(i);
    prior_precision(:,:,j) = inv(prior_precision(:,:,j));
    prior_mean(:,j)        = prior_precision(:,:,j) * prior_mean(:,j);
end

for i=1:ID(end)
    index   = same_indv{i};
    x_var_x = zero_mat_x;
    x_var_y = zero_mat_y;
    x_var_x = x_var_x + xr(index,:)' * (precision(i)) * xr(index,:);
    x_var_y = x_var_y + xr(index,:)' * (precision(i)) * ...
              (y_minus_xf_alpha(index) - r_int(i));
        
    post_cov_mat = inv(x_var_x + prior_precision(:,:, grp(i)));
    post_mean    = post_cov_mat * (x_var_y + prior_mean(:, grp(i)));
    cholesky     = chol(post_cov_mat, 'lower');
    output(:,i)  = cholesky * random(i,:)' + post_mean;
end
end