%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Article : Multiclass classification of growth curves using random change 
%            points and heterogeneous random effects
%
%  Authors : Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and
%            Scott A. Sisson
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import the data
data = readmatrix("../data/rotavirus.csv");

% Indicator for variable knot
variable_knot = 1;

% Number of knots
nknot = 3;

% Number of MCMC iterations
n_mcmc = 200000;

% Number of burnin
burnin = 100000;

% Thinning
thin   = 25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prior distribution parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_pred_xf      = 1;
prior_alpha.mat  = 25 * eye(num_pred_xf);
prior_alpha.mean = zeros(num_pred_xf, 1);

variance_scale  = 5;

prior_var.delta = 1;
prior_var.A     = variance_scale;
prior_var.a     = 1;
prior_var.Xi    = 1;

prior_rint.delta = 1;
prior_rint.scale = variance_scale;

num_pred_xr        = 4;
base_dist.mean_vec = zeros(num_pred_xr, 1);
base_dist.scale    = 1 * eye(num_pred_xr);
base_dist.lambda   = 0.001;
base_dist.df       = num_pred_xr + 1;

prior_eta.v     = 2;
prior_eta.delta = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result   = cell(1, 1);
init_grp = 2;

time = (1:nknot)/(nknot+1);

rng(1132029);
init.weight      = ones(init_grp,1) / init_grp;
init.grp_cov_mat = zeros(num_pred_xr, num_pred_xr, init_grp);

for i=1:init_grp
    init.grp_cov_mat(:,:,i) = iwishrnd(base_dist.scale, base_dist.df);
end

init.grp_mean = mvnrnd(base_dist.mean_vec, init.grp_cov_mat /...
                       base_dist.lambda, init_grp)';
init.eta      = 1;
init.knot     = time;

grp    = zeros(373,1);
beta_i = zeros(num_pred_xr, 373);

for i=1:373
    grp(i)      = randsample(init_grp, 1, true, init.weight);
    beta_i(:,i) = mvnrnd(init.grp_mean(:,grp(i)), ...
                  init.grp_cov_mat(:,:,grp(i)));
end

init.grp    = grp;
init.beta_i = beta_i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparing data and running MCMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obs_time = data(:,2);
ID       = data(:,1);
y        = data(:,3);
xf       = ones(length(ID),1);
xr       = obs_time;
        
for i=1:nknot
    tmp = obs_time - time(i);
    negative = tmp < 0;
    tmp(negative) = 0;
    xr = [xr, tmp];
end
for i=1:(size(xr,2)-1)
    xr(:,i) = xr(:,i) - xr(:,(i+1));
end

same_indv = cell(ID(end),1);
for i=1:ID(end)
    same_indv{i} = find(ID == i);
end
Ti = cellfun('length', same_indv);
        
[result{1}.alpha, result{1}.sigma2, result{1}.beta_i, result{1}.r_int, ...
    result{1}.var_reff, result{1}.grp_mean, result{1}.grp_cov_mat, ...
    result{1}.knot, result{1}.grp, result{1}.active_grp, ...
    result{1}.weight, result{1}.eta] = mcmc(n_mcmc, ID, y, xf, xr, ...
        obs_time, same_indv, Ti, base_dist, prior_alpha, prior_var, ...
        prior_rint, prior_eta, init, variable_knot);

est_grp_label = result{1}.grp(:,((burnin+1):thin:n_mcmc));

writematrix(est_grp_label, "../results/application/est_grp_label.csv");

est_alpha = result{1}.alpha(:,((burnin+1):thin:n_mcmc));
est_beta_i = result{1}.beta_i(:, :, ((burnin+1):thin:n_mcmc));
est_r_int = result{1}.r_int(:,((burnin+1):thin:n_mcmc));
est_knot = result{1}.knot(:, :, ((burnin+1):thin:n_mcmc));
save("../results/application/est_alpha.mat", "est_alpha");
save("../results/application/est_beta_i.mat", "est_beta_i");
save("../results/application/est_r_int.mat", "est_r_int");
save("../results/application/est_knot.mat", "est_knot");
