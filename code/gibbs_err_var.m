%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Article : Multiclass classification of growth curves using random change 
%            points and heterogeneous random effects
%
%  Authors : Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and
%            Scott A. Sisson
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gibbs sampling for variance of errors
function [sigma2, omega] = gibbs_err_var(ID, y , xf_alpha, xr, beta_i, ...
                                         r_int, old_sigma2, old_omega, ...
                                         same_indv, Ti, prior_var)

delta            = prior_var.delta;
A                = prior_var.A;
a                = prior_var.a;
Xi               = prior_var.Xi;
error_sq         = zeros(ID(end), 1);
y_minus_xf_alpha = y - xf_alpha;

for i=1:ID(end)
    index       = same_indv{i};
    error_sq(i) = sum((y_minus_xf_alpha(index) - r_int(i) - ...
                  xr(index,:) * beta_i(:,i)).^2);
end

a_sigma2 = 1 ./ gamrnd(0.5, 1 / (A^(-2) + delta / old_sigma2));
sigma2   = 1 ./ gamrnd((delta + length(ID)) / 2, 1 / (delta / a_sigma2 +...
           1 / 2 * sum(error_sq .* old_omega)));

omega = gamrnd((Ti / 2 + a), 1 ./ (Xi + 1 / (2 * sigma2) * error_sq));
end