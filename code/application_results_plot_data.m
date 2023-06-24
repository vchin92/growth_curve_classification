%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Article : Multiclass classification of growth curves using random change 
%            points and heterogeneous random effects
%
%  Authors : Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and
%            Scott A. Sisson
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import data
data = readmatrix("../data/rotavirus.csv");

% Load estimation results from MCMC
r_int = load(strcat('../results/application/est_r_int.mat'), 'est_r_int');
r_int = r_int.('est_r_int');
r_int = mean(r_int,2);

alpha = load(strcat('../results/application/est_alpha.mat'), 'est_alpha');
alpha = alpha.('est_alpha');
alpha = mean(alpha);

knot = load(strcat('../results/application/est_knot.mat'), 'est_knot');
knot = knot.('est_knot');

beta_i = load(strcat('../results/application/est_beta_i.mat'), 'est_beta_i');
beta_i = beta_i.('est_beta_i');

% Store knot locations and calculate slopes before and after knot changes
[P , n_dim , n_sample] = size(knot);
ordered_param = cell(n_dim,P);

for i=1:P
    for j=1:n_dim
        unique_knot = sort(unique(knot(i,j,:)));
        segment_param = zeros(n_dim-1, length(unique_knot)+1);

        for k=1:(length(unique_knot)+1)
            
            if k==1
                slope = mean([squeeze(beta_i(j,i, knot(i,j,:) >= ...
                              unique_knot(k))); squeeze(beta_i(j+1,i, ...
                              knot(i,j,:) < unique_knot(k)))]);
            else
                slope = mean([squeeze(beta_i(j,i, knot(i,j,:) > ...
                              unique_knot(k-1))); squeeze(beta_i(j+1,i, ...
                              knot(i,j,:) <= unique_knot(k-1)))]);
            end

            if k==(length(unique_knot)+1)
                knot_location = j/3;
            else
                knot_location = unique_knot(k);
            end

            segment_param(:,k) = [knot_location, slope];
        end
        ordered_param{j,i} = segment_param;
    end
end

rgb = [51 0 0;
    0 0 255;
    255 0 0;
    128 128 128;
    0 255 0;
    0 255 255;
    255 0 255]/255;

% Load estimated clusters
est_cluster = readmatrix("../results/application/cluster.csv");

cluster = cell(1,max(est_cluster));
for i=1:max(est_cluster)
    cluster{1,i} = find(est_cluster==i);
end

same_indv = cell(P,1);
for i=1:P
    same_indv{i} = find(data(:,1) == i);
end

% Plotting Figure 4
rng(1132029);

figure(1)
for i=1:length(cluster)
    subplot(2,4,i)
    tmp = cluster{1,i};

    if length(tmp) > size(rgb,1)
        selected{i} = sort(randsample([cluster{1,i}],size(rgb,1)));
    else
        selected{i} = cluster{1,i};
    end

    index = same_indv{selected{i}(1)};   
    plot(data(index,2), data(index,3),'linewidth',1.5,'color',rgb(1,:))
    ylim([-6 6])
    box off
    hold on

    for j=2:(min(length(cluster{1,i}),size(rgb,1)))
        index = same_indv{selected{i}(j)};
        plot(data(index,2), data(index,3),'linewidth',1.5,'color',rgb(j,:))
    end

    title(strcat('Subgroup', {' '}, num2str(i),' (', ...
          num2str(length(cluster{i})), ')'))
    xlabel('$t$','Interpreter','latex')
    ylabel('HAZ')
    set(gca,'FontSize',20)
end

cd ../figures
saveas(gcf,'figure4','epsc')

% Plotting Figure 5
figure(2)
for i=1:length(cluster)
    subplot(2,4,i)
    knot_dist = cell(1,n_dim);

    for k=1:n_dim
        knot_dist{1,k} = diff([(k-1)/3 ...
                               ordered_param{k,selected{i}(1)}(1,:)]);
    end
    
    z = [r_int(selected{i}(1)) + alpha, ...
            knot_dist{1,1} .* ordered_param{1,selected{i}(1)}(2,:), ...
            knot_dist{1,2} .* ordered_param{2,selected{i}(1)}(2,:), ...
            knot_dist{1,3} .* ordered_param{3,selected{i}(1)}(2,:)];
     
    x  = [0 ordered_param{1,selected{i}(1)}(1,:) ...
            ordered_param{2,selected{i}(1)}(1,:) ...
            ordered_param{3,selected{i}(1)}(1,:)]';
    
    plot(x, cumsum(z),'linewidth',1.5,'color',rgb(1,:))
    ylim([-6 6])
    box off
    hold on
    for j=2:length(selected{i})
        knot_dist = cell(1,n_dim);

        for k=1:n_dim
            knot_dist{1,k} = diff([(k-1)/3 ...
                                    ordered_param{k,selected{i}(j)}(1,:)]);
        end
        
        z = [r_int(selected{i}(j)) + alpha, ...
                knot_dist{1,1} .* ordered_param{1,selected{i}(j)}(2,:), ...
                knot_dist{1,2} .* ordered_param{2,selected{i}(j)}(2,:), ...
                knot_dist{1,3} .* ordered_param{3,selected{i}(j)}(2,:)];
       
        x  = [0 ordered_param{1,selected{i}(j)}(1,:) ...
                ordered_param{2,selected{i}(j)}(1,:) ...
                ordered_param{3,selected{i}(j)}(1,:)]';

        plot(x, cumsum(z),'linewidth',1.5,'color',rgb(j,:))
    end

    title(strcat('Subgroup', {' '}, num2str(i),' (', ...
            num2str(length(cluster{i})), ')'))
    xlabel('$t$','Interpreter','latex')
    ylabel('HAZ')
    set(gca,'FontSize',20)
end

saveas(gcf,'figure5','epsc')