%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Article : Multiclass classification of growth curves using random change 
%            points and heterogeneous random effects
%
%  Authors : Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and
%            Scott A. Sisson
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gibbs sampling for group labels
function [grp_no, n_grp] = gibbs_label(ID, y, xf_alpha, xr, r_int, ...
                                       grp_mean, grp_cov_mat, precision,...
                                       u, weight, same_indv)

y_minus_xf_alpha = y - xf_alpha;
neg_inf_vec      = -Inf * ones(size(grp_mean, 2), 1);
grp_no           = zeros(ID(end), 1);
n_grp            = zeros(size(grp_mean, 2), 1);

for i=1:ID(end)
    index = same_indv{i};
    Ti    = length(index);
    Z     = y_minus_xf_alpha(index) - r_int(i);
    XR    = xr(index,:);
    
    log_lik = neg_inf_vec;
    class   = find(weight > u(i));
    
    if length(class) > 1
        for j=1:length(class)
            log_lik(class(j)) = sum(log(normpdf(Z, XR * grp_mean(:, ...
                                class(j)), sqrt(diag(XR * grp_cov_mat(:,...
                                :, class(j)) * XR' + 1 / precision(i) * ...
                                eye(Ti))))));      
        end
        grp_no(i)        = randsample(class, 1, true, exp(log_lik(class)));
        n_grp(grp_no(i)) = n_grp(grp_no(i)) + 1;
    else
        grp_no(i)    = class;
        n_grp(class) = n_grp(class) + 1;
    end
end
end