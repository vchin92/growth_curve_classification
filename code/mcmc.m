%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Article : Multiclass classification of growth curves using random change 
%            points and heterogeneous random effects
%
%  Authors : Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and
%            Scott A. Sisson
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main function for running MCMC
function [alpha, sigma2, beta_i, r_int, var_reff, grp_mean, grp_cov_mat,...
    knot, grp, active_grp, weight, eta] = mcmc(n_mcmc, ID, y, xf, xr, ...
    obs_time, same_indv, Ti, base_dist, prior_alpha, prior_var, ...
    prior_rint, prior_eta, init, variable_knot)

rng(1132029);

num_pred_xf = size(xf, 2);

alpha  = zeros(num_pred_xf, n_mcmc);
sigma2 = ones(1, n_mcmc);
omega  = ones(ID(end), 1);
eta    = init.eta * ones(1, n_mcmc);

num_pred_xr   = size(xr, 2);
beta_i        = zeros(num_pred_xr, ID(end), n_mcmc);
beta_i(:,:,1) = init.beta_i;

grp      = zeros(ID(end), n_mcmc);
grp(:,1) = init.grp;

active_grp    = zeros(1, n_mcmc);
active_grp(1) = length(unique(grp(:,1)));

init_grp                 = length(init.weight);
total                    = 60;
grp_mean                 = zeros(num_pred_xr, total, n_mcmc);
grp_mean(:,1:init_grp,1) = init.grp_mean;

grp_cov_mat                   = zeros(num_pred_xr, num_pred_xr, total, ...
                                n_mcmc);
grp_cov_mat(:,:,1:init_grp,1) = init.grp_cov_mat;

weight               = zeros(total, n_mcmc);
weight(1:init_grp,1) = init.weight;

r_int    = zeros(ID(end), n_mcmc);
var_reff = ones(1, n_mcmc);

knot        = zeros(ID(end), num_pred_xr-1, n_mcmc);
knot(:,:,1) = init.knot .* ones(ID(end), num_pred_xr-1);

for i=2:n_mcmc
    
    % Mean of random intercepts
    precision  = omega / sigma2(i-1);
    alpha(:,i) = gibbs_alpha(ID, y, xf, xr, beta_i(:,:,i-1), ...
                             r_int(:,i-1), precision, same_indv, ...
                             prior_alpha);
    
    % Variance of errors
    xf_alpha  = xf * alpha(:,i);
    sigma2(i) = gibbs_err_var(ID, y , xf_alpha, xr, beta_i(:,:,i-1), ...
                              r_int(:,i-1), sigma2(i-1), omega, ...
                              same_indv, Ti, prior_var);
    
    % Random slopes
    precision     = omega / sigma2(i);
    beta_i(:,:,i) = gibbs_beta_i(ID, y, xf_alpha, xr, r_int(:,i-1), ...
                                 precision, grp(:,i-1), ...
                                 grp_mean(:,:,i-1), ...
                                 grp_cov_mat(:,:,:,i-1), same_indv);
    
    % Random intercept and its variance
    [r_int(:,i), var_reff(i)] = gibbs_rint(ID, y , xf_alpha, xr, ...
                                           beta_i(:,:,i), precision, ...
                                           var_reff(i-1), same_indv, Ti,...
                                           prior_rint);
                                 
    % Group parameters
    [m, c] = gibbs_grp_param(grp(:,i-1), grp_mean(:,:,i-1), ...
                             grp_cov_mat(:,:,:,i-1), beta_i(:,:,i), ...
                             base_dist);
                         
    grp_mean(:,:,i)      = m;
    grp_cov_mat(:,:,:,i) = c;
    
    % Auxiliary variable for random truncation
    u = gibbs_u(grp(:,i-1), weight(:,i-1));
    
    if i > 2
        k            = find(weight(:,i-1)>0, 1, 'last' );
        trunc_weight = weight(k, i-1);
        u_min        = min(u);
        w            = weight(1:k, i-1);
        while u_min < trunc_weight
            stick = betarnd(1, eta(i-1));
            w     = [w(1:(k-1)); trunc_weight * stick; trunc_weight * ...
                    (1 - stick)];
            k     = k + 1;
            
            new_cov              = iwishrnd(base_dist.scale, base_dist.df);
            grp_cov_mat(:,:,k,i) = new_cov;
            grp_mean(:,k,i)      = mvnrnd(base_dist.mean_vec, new_cov / ...
                                   base_dist.lambda);

            trunc_weight = w(end);
        end
        weight(1:length(w),i-1) = w;
    end
    
    % Random knots
    if variable_knot == 1
        [xr, knot(:,:,i)] = mh_knot(y, xf_alpha, xr, r_int(:,i), ...
                                    beta_i(:,:,i), precision, ...
                                    knot(:,:,i-1), same_indv, Ti, obs_time);
    end
    
    % Group label
    [grp(:,i), n_grp] = gibbs_label(ID, y, xf_alpha, xr, r_int(:,i), ...
                                    grp_mean(:,:,i), ...
                                    grp_cov_mat(:,:,:,i), precision, u, ...
                                    weight(:,i-1), same_indv);
    active_grp(i) = length(unique(grp(:,i)));                             
    
    % Rearranging group parameters
    non_empty   = unique(grp(:,i));
    first_empty = find(n_grp < 1, 1);
    if first_empty < max(non_empty)
        for k=first_empty:length(non_empty)
            re_index        = grp(:,i) == non_empty(k);
            grp(re_index,i) = k;
        end
    end
       
    n_grp                                = n_grp(non_empty);
    tmp_g_mean                           = grp_mean(:,non_empty,i);
    tmp_g_cov                            = grp_cov_mat(:,:,non_empty,i);
    tmp_weight                           = weight(non_empty,i-1);
    weight(1:length(non_empty),i-1)      = tmp_weight;
    weight(length(non_empty)+1:end, i-1) = 0;
                            
    % Weights
    [tmp_w, tmp_m, tmp_c] = gibbs_weight(n_grp, tmp_g_mean, tmp_g_cov, ...
                                         u, eta(i-1), base_dist);
    
    l_w = length(tmp_w);
    l_e = total - l_w;
    
    weight(1:l_w,i)              = tmp_w;
    grp_mean(:,1:l_w,i)          = tmp_m;
    grp_mean(:,l_w+1:end,i)      = zeros(num_pred_xr, l_e);
    grp_cov_mat(:,:,1:l_w,i)     = tmp_c;
    grp_cov_mat(:,:,l_w+1:end,i) = repmat(zeros(num_pred_xr, ...
                                   num_pred_xr), [1 1 l_e]);
                               
    % Concentration parameter
    eta(i) = gibbs_eta(grp(:,i), eta(i-1), prior_eta);

    disp(i)
end                               
end