%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Article : Multiclass classification of growth curves using random change 
%            points and heterogeneous random effects
%
%  Authors : Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and
%            Scott A. Sisson
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gibbs sampling for random intercept and variance of normal distribution
% describing random intercept
function [r_int, var_reff] = gibbs_rint(ID, y , xf_alpha, xr, beta_i, ...
                                        precision, old_var_reff, ...
                                        same_indv, Ti, prior_rint)

delta            = prior_rint.delta;
A                = prior_rint.scale;
error            = zeros(ID(end), 1);
y_minus_xf_alpha = y - xf_alpha;

for i=1:ID(end)
    index    = same_indv{i};
    error(i) = sum(y_minus_xf_alpha(index) - xr(index,:) * beta_i(:,i));
end

post_variance = 1 ./ (Ti .* precision + 1 / old_var_reff);
post_mean     = post_variance .* precision .* error;
r_int         = randn(ID(end), 1) .* sqrt(post_variance) + post_mean;

a_sigma_reff = 1 ./ gamrnd(0.5, 1 / (A^(-2) + delta / old_var_reff));
var_reff     = 1 ./ gamrnd((delta + ID(end)) / 2, 1 / (delta / ...
               a_sigma_reff + 1 / 2 * sum(r_int .^2)));
end