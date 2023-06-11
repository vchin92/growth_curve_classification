%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Article : Multiclass classification of growth curves using random change 
%            points and heterogeneous random effects
%
%  Authors : Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and
%            Scott A. Sisson
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Metropolis-Hastings algorithm for sampling random knots
function [old_xr, old_knot] = mh_knot(y, xf_alpha, old_xr, r_int, ...
                                      beta_i, precision, old_knot, ...
                                      same_indv, Ti, obs_time)

nknot    = size(old_knot, 2);
zero_vec = zeros(size(old_knot, 1), 1);
one_vec  = ones(size(old_knot, 1), 1);
n_ppl    = length(same_indv);
random   = rand(1, n_ppl * nknot);
m        = 1;
old_ll   = zero_vec;
Z        = y - xf_alpha - repelem(r_int, Ti);
sim_knot = rand(size(beta_i, 2), nknot) * 1/nknot + [0:(nknot-1)]/nknot;

for j=1:nknot
    new_knot = old_knot;
    new_xr   = old_xr;
    
    tmp_old_knot  = [zero_vec, old_knot, one_vec];
    new_knot(:,j) = sim_knot(:,j); 
    tmp_new_knot  = [zero_vec, new_knot, one_vec];
    
    d_new = diff(tmp_new_knot, 1, 2);
    d_old = diff(tmp_old_knot, 1, 2);
    r_new = d_new(:,j+1) .* d_new(:,j);
    r_old = d_old(:,j+1) .* d_old(:,j);
        
    if j==1
        tmp1 = obs_time - positive_part(obs_time - repelem(new_knot(:,j), Ti));
    else
        tmp1 = positive_part(obs_time - repelem(new_knot(:,j-1), Ti)) - ...
               positive_part(obs_time - repelem(new_knot(:,j), Ti));
    end
    
    if j~=nknot
        tmp2 = positive_part(obs_time - repelem(new_knot(:,j), Ti)) - ...
               positive_part(obs_time - repelem(new_knot(:,j+1), Ti));
    else
        tmp2 = positive_part(obs_time - repelem(new_knot(:,j), Ti));
    end
    
    new_xr(:,j:(j+1)) = [tmp1, tmp2];
    
    for i=1:size(new_knot, 1)
        index  = same_indv{i};
        NEW_XR = new_xr(index,:);
        new_ll = sum(log(normpdf(Z(index) - NEW_XR * beta_i(:,i), 0, ...
                 sqrt(1 / precision(i)))));
        if j==1
            XR        = old_xr(index,:);
            old_ll(i) = sum(log(normpdf(Z(index) - XR * beta_i(:,i), 0, ...
                        sqrt(1 / precision(i)))));
        end
        
        if exp(new_ll - old_ll(i)) * r_new(i) / r_old(i) > random(m)
            old_ll(i) = new_ll;
            old_knot(i,j) = new_knot(i,j);
            old_xr(index,:) = new_xr(index,:);
        end
        m = m + 1;
    end
end
end