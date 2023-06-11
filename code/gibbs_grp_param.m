%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Article : Multiclass classification of growth curves using random change 
%            points and heterogeneous random effects
%
%  Authors : Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and
%            Scott A. Sisson
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gibbs sampling for mean vectors and covariance matrices of normal
% distribution describing random slopes for each subgroup
function [grp_m, grp_c] = gibbs_grp_param(grp, grp_mean, grp_cov_mat, ...
                                          beta_i, base_dist)

active_grp = unique(grp);
same_grp   = cell(1, length(active_grp));
grp_c      = grp_cov_mat;
grp_m      = grp_mean;
mean_vec   = base_dist.mean_vec;
lambda     = base_dist.lambda;
df         = base_dist.df;
scale      = base_dist.scale;

for i=1:length(active_grp)
    j            = active_grp(i);
    same_cluster = find(grp == j);
    same_grp{i}  = same_cluster;
    n_grp        = length(same_cluster);
    
    mean_sample = mean(beta_i(:, same_cluster), 2);
    mu_n        = 1 / (lambda + n_grp) * (lambda * mean_vec + n_grp * ...
                  mean_sample);
    l_n         = lambda + n_grp;
    v_n         = df + n_grp;
    difference  = (beta_i(:, same_cluster) - mean_sample);
    cov_mat     = difference * difference';
    psi_n       = scale + cov_mat + n_grp * lambda / (n_grp + lambda) * ...
                  ((mean_sample - mean_vec) * (mean_sample - mean_vec)');
                
    grp_c(:,:,j) = iwishrnd(psi_n, v_n);
    grp_m(:,j)   = mvnrnd(mu_n, grp_c(:,:,j) / l_n);
end
end