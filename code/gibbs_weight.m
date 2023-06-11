%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Article : Multiclass classification of growth curves using random change 
%            points and heterogeneous random effects
%
%  Authors : Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and
%            Scott A. Sisson
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gibbs sampling for weights
function [w, grp_m, grp_cov] = gibbs_weight(n_grp, grp_mean, ...
                                            grp_cov_mat, u, eta, base_dist)

mean_vec   = base_dist.mean_vec;
lambda     = base_dist.lambda;
df         = base_dist.df;
scale      = base_dist.scale;

active_grp    = [n_grp > 0; true];
w             = zeros(length(active_grp), 1);
random_w      = gamrnd([n_grp(active_grp(1:end-1)); eta], 1);
random_w      = random_w / sum(random_w);
w(active_grp) = random_w;

grp_cov        = grp_cov_mat;
grp_m          = grp_mean;
k              = length(w);
new_cov        = iwishrnd(scale, df);
grp_cov(:,:,k) = new_cov;
grp_m(:,k)     = mvnrnd(mean_vec, new_cov / lambda);

u_min        = min(u);
trunc_weight = w(end);

while u_min < trunc_weight
    stick = betarnd(1, eta);
    w     = [w(1:(k-1)); trunc_weight * stick; trunc_weight * (1 - stick)];
    k     = k + 1;
    
    new_cov        = iwishrnd(scale, df);
    grp_cov(:,:,k) = new_cov;
    grp_m(:,k)     = mvnrnd(mean_vec, new_cov / lambda);
    
    trunc_weight = w(end);
end  
end