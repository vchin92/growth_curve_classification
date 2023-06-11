%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Article : Multiclass classification of growth curves using random change 
%            points and heterogeneous random effects
%
%  Authors : Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and
%            Scott A. Sisson
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script for generating simulated data and plotting Figure 1

df_fixed    = cell(100, 8); % Dataset with fixed knot locations
df_random   = cell(100, 8); % Dataset with random knot locations
knot_random = cell(100, 8); % Knot locations for df_random
label       = cell(100, 8); % True group label
beta        = cell(100, 8); % Random slopes
r_eff       = cell(100, 8); % Random intercepts

for nsize=1:8
    for iter=1:100

        rng(1132000+iter)

        P         = 50+(nsize-1)*50; % Number of individuals
        weight    = [0.25; 0.25; 0.25; 0.25]; % Weight of mixture dist.
        sd_error  = sqrt(0.15); % Std deviation of Gaussian random errors
        t_min     = 10; % Min number of observations
        t_max     = 20; % Max number of observations
        reff_mean = 0.75; % Mean of Gaussian random intercepts
        reff_var  = 0.5; % Variance of Gaussian random intercepts
        grp_mean  = [-3, -7.5, -3,  4;
                     -3, -5, -1,  1;
                     -3,  0,  3, -3]; % Mean of Gaussian random slopes
        scale_mat = [1, -0.5, 0;
                     -0.5, 1, -0.5; 
                     0, -0.5, 1]; % Scale matrix of inverse-Wishart dist.


        index = randsample(1:4, P, true, weight);
        P_grp = tabulate(index);
        P_grp = P_grp(:,2); % Number of individuals in each group
        Ti    = randsample(t_min:t_max, P, true)';
        ID    = repelem(1:P, 1, Ti)'; % Subject ID for each individual


        label{iter, nsize} = sort(index);
        reff               = randn(ID(end),1)*sqrt(reff_var)+reff_mean;
        r_eff{iter, nsize} = reff;    
        error              = randn(length(ID),1)*sd_error; % Random errors

        
        % Generate random slopes
        P_i         = [0; cumsum(P_grp)];
        beta_i      = zeros(size(grp_mean, 1), P);
        grp_cov_mat = repmat(0.2*corrcov(iwishrnd(scale_mat, 100)), ...
                      [1 1 4]);
        for i=1:4
            beta_i(:,(P_i(i)+1):P_i(i+1)) = mvnrnd(grp_mean(:,i), ...
                                        grp_cov_mat(:,:,i), P_grp(i))';
        end
        beta{iter, nsize} = beta_i;

        
        % Generate dataset
        same_indv = cell(ID(end),1);
        for i=1:ID(end)
            same_indv{i} = find(ID==i);
        end

        xr_tmp = [];
        for i=1:P
            xr_tmp = [xr_tmp; sort(randsample(1:365, Ti(i), false))'];
        end

        xr_tmp   = xr_tmp/365;
        obs_time = xr_tmp;

        for knot_val=1:2
            xr          = xr_tmp;
            random_knot = knot_val;
            knot        = zeros(P, 2);
            if random_knot==1
                for i=1:P
                    knot(i,:) = [rand(1)*0.5, rand(1)*0.5+0.5]; % Random knot
                end
            else
                knot = repelem([1/3, 2/3], P, 1); % Fixed knot
            end

            for i=1:2
                tmp           = obs_time-repelem(knot(:,i), Ti, 1);
                negative      = tmp<0;
                tmp(negative) = 0;
                xr            = [xr, tmp];
            end

            for i=1:(size(xr,2)-1)
                xr(:,i) = xr(:,i) - xr(:,(i+1));
            end

            y = [];
            for i=1:P
                y = [y; reff(i)+xr(same_indv{i},:)*beta_i(:,i)+ ...
                    error(same_indv{i})];
            end
            data = [ID, obs_time, y];

            if random_knot == 1
                df_random{iter, nsize} = data;
                knot_random{iter, nsize} = knot;
            else
                df_fixed{iter, nsize} = data;
            end
        end
    end
end

cd ../data
save("df_fixed.mat", "df_fixed");
save("df_random.mat", "df_random");

% Plotting Figure 1
rng(1132000);

sample = [randsample(find(label{iter,nsize}==1),1);
          randsample(find(label{iter,nsize}==2),1);
          randsample(find(label{iter,nsize}==3),1);
          randsample(find(label{iter,nsize}==4),1)];

figure(1)
% Samples from dataset with fixed knot locations
for i=1:4
    subplot(2,4,i)
    plot(df_fixed{iter,nsize}(same_indv{sample(i)},2), ...
        df_fixed{iter,nsize}(same_indv{sample(i)},3), ...
        'x','linewidth',1.5,'Color','black')
    xlim([0 1])
    ylim([-6 6])
    xlabel('$t$','Interpreter','latex')
    ylabel('HAZ')
    set(gca,'FontSize',20)
    title(strcat('Subgroup', {' '}, num2str(i)))
    box off
    hold on
    plot([0 1/3 2/3 1], cumsum([r_eff{iter,nsize}(sample(i)) ...
        beta{iter,nsize}(:,sample(i))'*1/3]),'linewidth',1.5,'Color','black')
    plot([1/3 1/3], [-6 6],'--','Color','black')
    plot([2/3 2/3], [-6 6],'--','Color','black') 
end

% Samples from dataset with random knot locations
for i=1:4
    subplot(2,4,i+4)
    plot(df_random{iter,nsize}(same_indv{sample(i)},2), ...
        df_random{iter,nsize}(same_indv{sample(i)},3), ...
        'x','linewidth',1.5,'Color','black')
    xlim([0 1])
    ylim([-6 6])
    box off
    xlabel('$t$','Interpreter','latex')
    ylabel('HAZ')
    set(gca,'FontSize',20)
    title(strcat('Subgroup', {' '}, num2str(i)))
    hold on
    plot([0 knot_random{iter,nsize}(sample(i),:) 1], ...
        cumsum([r_eff{iter,nsize}(sample(i)) ...
        beta{iter,nsize}(:,sample(i))' .* ...
        diff([0 knot_random{iter,nsize}(sample(i),:) 1])]), ...
        'linewidth',1.5,'Color','black')
    plot([knot_random{iter,nsize}(sample(i),1) ...
        knot_random{iter,nsize}(sample(i),1)], [-6 6],'--','Color','black')
    plot([knot_random{iter,nsize}(sample(i),2) ...
        knot_random{iter,nsize}(sample(i),2)], [-6 6],'--','Color','black')
end

cd ../figures
saveas(gcf,'figure1','epsc')