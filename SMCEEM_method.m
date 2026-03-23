function [param_vals, prior_sample] = SMCEEM_method(args,prior_funcs)

% The Sequential Monte Carlo - Ensemble Ecosystem Modelling method proposed 
% in Vollert, Drovandi & Adams (2024). Unlocking ensemble ecosystem
% modelling for large and complex networks. PLOS Computational Biology
% 20(3) e1011976. https://doi.org/10.1371/journal.pcbi.1011976
%  
% Parameters: 
%   prior_funcs- data structure for the prior distribution of the parameters. 
%               Requires the following fields. 
%               sampler    - a sampling function that generates random vectors from
%                            the joint prior distribution;
%               pdf        - the joint probability density function;
%               trans_f    - transform of prior paramete space to ensure 
%                            unbounded support for MCMC sampling.
%               trans_finv - inverse of transform.
%  args       - data structure for arguments
%               skip_parameters - an array defining which parameters in the
%                                 interaction matrix are zero 
%               n_species       - the number of species in the system
%               num_params      - the number of parameters in the system
%               lower           - the prior lower bounds of the parameters
%               upper           - the prior upper bounds of the parameters.
%               n_particles     - number of particles for SMC sampler.
% Returns:
%    param_vals     - the parameter values for the EEM
%    prior_sample   - the prior sample of parameter values
%

%% Specify tuning parameters
a = 0.6;                % adaptive selection of discrepancy threshold sequence
c = 0.01;               % the goal for percentage unmoved particles in MCMC.
p_acc_min = 0.0001;     % minimum  acceptance rate in the MCMC interations before stopping.  
mcmc_trials = 10;       % the number of MCMC steps to try before selecting appropriate number
dist_final = 0;         % the target discrepancy threshold.

%% Begin Sampling

% values for adaptive steps
num_drop = floor(args.n_particles*a);                                            %Number of particles dropped each iteration
num_keep = args.n_particles-num_drop;                                            %Number of particles kept each iteration

% Initialise particle data structures
param_vals = zeros(args.n_particles,args.num_params);                            %parameter values for each particle
param_disc = Inf*ones(args.n_particles,1);                                          %discrepancy metric for each particle

% initial prior rejection algorithm
while sum(param_disc==Inf)>0
    parfor i = 1:args.n_particles
        if param_disc(i) ==Inf
        % sample prior
        param_vals(i,:) = prior_funcs.sampler(args);                            %part vals is not transformed
        % simulate model
        param_disc(i) = args.simulate_and_calc_discrepancy(param_vals(i,:),args);
        end
    end
end

%save prior sample
prior_sample = param_vals;

% transform the parameters
for i=1:args.n_particles
    param_vals(i,:) = prior_funcs.trans_f(param_vals(i,:),args);            %part vals is transformed
end

% sort the particles
[param_disc,ix] = sort(param_disc); 
param_vals = param_vals(ix,:); 

% determine next disprepacy threshold
dist_max = param_disc(args.n_particles);
dist_next = param_disc(num_keep);
fprintf('Current maximum discrepancy %f, now trying for %f, want to get to %f\n',dist_max,dist_next,dist_final)

% interate toward target discrepancy
while (dist_max > dist_final)

    %RESAMPLE
    % compute the covariance matrix (of particles that remain)
    cov_matrix = (2.38^2)*cov(param_vals(1:num_keep,:))/size(param_vals,2);

    % resample
    r = randsample(num_keep, args.n_particles-num_keep, 'true');                 %duplicate good ones
    param_vals((num_keep+1):args.n_particles, :) = param_vals(r,:);
    param_disc((num_keep+1):args.n_particles) = param_disc(r);
    
    %track acceptances and simulations
    acc_counter = zeros(args.n_particles-num_keep,1);
    
    % MOVE
    parfor i = (num_keep+1):args.n_particles %for each particles to mov

        %run MCMC
        [param_vals(i,:), param_disc(i),acc_counter(i-num_keep)] = ...
            MCMC(mcmc_trials,param_vals(i,:),cov_matrix,args,...
            prior_funcs,param_disc(i),acc_counter(i-num_keep),dist_next)
    end

    % determine number of MCMC iterations to perfrom 
    acc_rate = sum(acc_counter)/(mcmc_trials*(args.n_particles-num_keep));       %calculate acceptance rate
    mcmc_iters =   floor(log(c)/log(1-acc_rate)+1);                         %predict how many mcmc steps are needed to get approx unchanged paricles c%
    
    % move step
    parfor i = (num_keep+1):args.n_particles
        %run MCMC
        [param_vals(i,:), param_disc(i),acc_counter(i-num_keep)] = ...
            MCMC(mcmc_iters-mcmc_trials,param_vals(i,:),cov_matrix,args,...
            prior_funcs,param_disc(i),acc_counter(i-num_keep),dist_next)
    end
 
    % calc percentage acceptance
    num_mcmc_iters = max(0, mcmc_iters - mcmc_trials) + mcmc_trials;
    p_acc = sum(acc_counter)/(num_mcmc_iters*(args.n_particles-num_keep));
    
    %update number of mcmc trials for next SMC step
    mcmc_trials = ceil(mcmc_iters/2);
        
    %%REWEIGHT
    % Compute the next distance and maximum distance
    [param_disc,ix] = sort(param_disc); param_vals = param_vals(ix,:); % sort the particles

    % if most of the particles are under the final target then don't drop
    if (sum((param_disc > dist_final)) < num_drop)
        num_drop = sum((param_disc > dist_final));
        num_keep = args.n_particles-num_drop;
    end
   
    %calculate distances
    dist_max = param_disc(args.n_particles);
    dist_next = param_disc(num_keep);
    
    % check to see if we will reach desired tolerance at next iteration
    if (dist_next < dist_final)
        dist_next = dist_final;
    end
    fprintf('The next distance is %f and the maximum distance is %f and the number to drop is %d\n',dist_next,dist_max,num_drop);

    % if acceptance is too low, then bail
    if (p_acc < p_acc_min)
        fprintf('Getting out as MCMC acceptance rate is below acceptable threshold\n');
        break;
    end

end %end SMC steps

%do the inverse transform
for i=1:args.n_particles
    param_vals(i,:) = prior_funcs.trans_finv(param_vals(i,:),args);
end

end



%% MCMC function
function [params, disc,i_acc] = MCMC(mcmc_steps,params,cov_matrix,args,prior_funcs,disc,i_acc,dist_next)
%This function runs the MCMC steps for a particle. 

for r = 1:mcmc_steps
    % Gaussian random walk
    part_vals_prop = mvnrnd(params,cov_matrix);

    % Transform back to calculate prior probs and discrepancy
    prop = prior_funcs.trans_finv(part_vals_prop, args);

    % Calculate prior probabilities
    prior_curr = prior_funcs.pdf(params,args);
    prior_prop = prior_funcs.pdf(part_vals_prop,args);

    % early rejection (assumes symmetric proposals)
    if (isnan(prior_prop/prior_curr) || rand > prior_prop/prior_curr)
        continue;
    end

    %find proposal discrepancy
    dist_prop = args.simulate_and_calc_discrepancy(prop,args);

    %Accept a particle if it is within the target distance.
    if (dist_prop <= dist_next)
        % then the metropolis-hastings ratio is satisfied
        params = part_vals_prop;
        disc = dist_prop;
        i_acc = i_acc + 1;
    end
end

end
