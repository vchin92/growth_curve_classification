%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Article : Multiclass classification of growth curves using random change 
%            points and heterogeneous random effects
%
%  Authors : Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and
%            Scott A. Sisson
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set data_type to 'df_fixed' for data generating model with fixed knots
% or 'df_random' for data generating model with random knots

% Set knot_type to 'fixed' for estimation model with fixed knots
% or 'random' for estimation model with random knots

data_type = 'df_fixed';
knot_type = 'fixed';

if strcmp(data_type, 'df_fixed')
    output_file_name1 = 'df';
else
    output_file_name1 = 'dr';
end

if strcmp(knot_type, 'fixed')
    variable_knot = 0;
    output_file_name2 = 'mf';
else
    variable_knot = 1;
    output_file_name2 = 'mr';
end


% Load data
data = load(strcat('../data/',data_type,'.mat'),data_type);
data = data.(data_type);

% Number of knots
nknot = 2;

% Number of MCMC iterations
n_mcmc = 100000;

% Number of burnin
burnin = n_mcmc/2;

% Thinning
thin   = 20;

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

num_pred_xr        = 3;
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
init     = cell(1, 8);
init_grp = 2;

time = (1:nknot)/(nknot+1);

for n_size=1:8
    rng(1132029);
    init{n_size}.weight      = ones(init_grp,1) / init_grp;
    init{n_size}.grp_cov_mat = zeros(num_pred_xr, num_pred_xr, init_grp);
    
    for i=1:init_grp
        init{n_size}.grp_cov_mat(:,:,i) = iwishrnd(base_dist.scale, ...
                                          base_dist.df);
    end
    
    init{n_size}.grp_mean = mvnrnd(base_dist.mean_vec, ...
                            init{n_size}.grp_cov_mat/base_dist.lambda, ...
                            init_grp)';
    init{n_size}.eta      = 1;
    init{n_size}.knot     = time;

    grp    = zeros(50+(n_size-1)*50,1);
    beta_i = zeros(num_pred_xr, 50+(n_size-1)*50);
    
    for i=1:(50+(n_size-1)*50)
        grp(i)      = randsample(init_grp, 1, true, init{n_size}.weight);
        beta_i(:,i) = mvnrnd(init{n_size}.grp_mean(:,grp(i)), ...
                      init{n_size}.grp_cov_mat(:,:,grp(i)));
    end
    
    init{n_size}.grp    = grp;
    init{n_size}.beta_i = beta_i;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparing data and running MCMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nsize=1:8

    for nrep = 1:100
    
    obs_time = data{nrep, nsize}(:,2);
    ID       = data{nrep, nsize}(:,1);
    y        = data{nrep, nsize}(:,3);
    xf       = ones(length(ID),1);
    xr       = data{nrep, nsize}(:,2);
            
    for i=1:nknot
        tmp = data{nrep, nsize}(:,2) - time(i);
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
            prior_rint, prior_eta, init{nsize}, variable_knot);
        
    est_grp_label = result{1}.grp(:,((burnin+1):thin:n_mcmc));
    
    writematrix(est_grp_label, strcat("../results/simulation/", ...
        output_file_name1,"_",output_file_name2,"_",num2str(nsize),"_", ...
        num2str(nrep),".csv"));  
    end
end