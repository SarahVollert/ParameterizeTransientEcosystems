function [param_vals, prior_sample] = SMCEEM_method_temporallyAdaptive(args,prior_funcs)

% This algorithm is a modification that incorporates temporal adaptations
% into the Sequential Monte Carlo - Ensemble Ecosystem Modelling method proposed 
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
b = 0.8;                % max particles to be kept when resampling
c = 0.01;               % the goal for percentage unmoved particles in MCMC.
d = 1;                % factor for increasing time step.
p_acc_min = 0.0001;     % minimum  acceptance rate in the MCMC interations before stopping.
mcmc_trials = 10;       % the number of MCMC steps to try before selecting appropriate number
dist_final = 0;         % the target discrepancy threshold.

T_sim0 = 0.01*args.T_max; %the initial period to simulate

%% Begin Sampling

%Time to simulate to
args.T_sim = T_sim0;

% Initialise particle data structures
param_vals = zeros(args.n_particles,args.num_params);                            %parameter values for each particle
param_disc = Inf*ones(args.n_particles,1);                                          %discrepancy metric for each particle

% initial prior rejection algorithm
while sum(param_disc==Inf)>0
    for i = 1:args.n_particles
        if param_disc(i) ==Inf

            % sample prior
            param_vals(i,:) = prior_funcs.sampler(args);                            %part vals is not transformed
            parameters = param_vals(i,:);

            %initialise simulation storage
            param_sims(i).x = 0; %saves simulation times
            param_sims(i).y = parameters(1:args.n_species)'; %saves simulation populations

            % simulate model and calculate discrepancy
            [param_disc(i), param_sims(i)] = args.simulate_and_calc_discrepancy(parameters,param_sims(i),[0 args.T_sim],args);
        end
    end
end
p_acc = 1; %probability of acceptance

%save prior sample
prior_sample = param_vals;

% transform the parameters
for i=1:args.n_particles
    param_vals_trans(i,:) = prior_funcs.trans_f(param_vals(i,:),args);            %part vals is transformed
end

% sort the particles
[param_disc,ix] = sort(param_disc);
param_vals_trans = param_vals_trans(ix,:);
param_sims = param_sims(ix);

% determine next disprepacy threshold
dist_max = param_disc(args.n_particles);
dist_next = param_disc((1-a)*args.n_particles);

% calculate number to keep
num_keep = max(sum(param_disc<=dist_next),floor(args.n_particles*(1-a)));   %Number of particles kept
num_drop = args.n_particles-num_keep;                                       %Number of particles dropped
fprintf('Simulation time=%f, max discrepancy=%f, next discrepancy=%f, number to drop=%d\n',args.T_sim,dist_max,dist_next,num_drop);

count_SMC_iters = 0;
%%
% interate toward target discrepancy
while (dist_max > dist_final || args.T_sim<args.T_max)

    count_SMC_iters = count_SMC_iters+1;
    %if we are dropping enough particles OR we have reached full simulation time
    if num_keep < floor(b*args.n_particles) || args.T_sim == args.T_max
        %RESAMPLE
        % compute the covariance matrix (of particles that remain)
        cov_matrix = (2.38^2)*cov(param_vals_trans(1:num_keep,:))/size(param_vals_trans,2);

        % resample
        r = randsample(num_keep, args.n_particles-num_keep, 'true');                 %duplicate good ones
        param_vals_trans((num_keep+1):args.n_particles, :) = param_vals_trans(r,:);
        param_disc((num_keep+1):args.n_particles) = param_disc(r);
        param_sims((num_keep+1):args.n_particles) = param_sims(r);

        
        %track acceptances and simulations
        acc_counter = zeros(args.n_particles-num_keep,1);

        % MOVE
        parfor i = (num_keep+1):args.n_particles %for each particles to mov

            %run MCMC
            [param_vals_trans(i,:), param_disc(i),acc_counter(i-num_keep),param_sims(i)] = ...
                MCMC(mcmc_trials,param_vals_trans(i,:),cov_matrix,args,...
                prior_funcs,param_disc(i),acc_counter(i-num_keep),dist_next,param_sims(i));
        end

        % determine number of MCMC iterations to perfrom
        acc_rate = sum(acc_counter)/(mcmc_trials*(args.n_particles-num_keep));       %calculate acceptance rate
        mcmc_iters =   floor(log(c)/log(1-acc_rate)+1);                         %predict how many mcmc steps are needed to get approx unchanged paricles c%

        % move step
        parfor i = (num_keep+1):args.n_particles
            %run MCMC
            [param_vals_trans(i,:), param_disc(i),acc_counter(i-num_keep),param_sims(i)] = ...
                MCMC(mcmc_iters-mcmc_trials,param_vals_trans(i,:),cov_matrix,args,...
                prior_funcs,param_disc(i),acc_counter(i-num_keep),dist_next,param_sims(i));
        end

        % calc percentage acceptance
        num_mcmc_iters = max(0, mcmc_iters - mcmc_trials) + mcmc_trials;
        p_acc = sum(acc_counter)/(num_mcmc_iters*(args.n_particles-num_keep));

        %update number of mcmc trials for next SMC step
        mcmc_trials = ceil(mcmc_iters/2);
    end

   
    %%REWEIGHT
    %increase the simulation time
    percent_nodisc = sum(param_disc==0)/args.n_particles;
    T_sim_old = args.T_sim;
    args.T_sim = min(T_sim_old*(1+d*percent_nodisc),args.T_max);

    %if we are increasing the simulation time
    if T_sim_old ~= args.T_sim
        parfor i=1:args.n_particles
            % transform bacl
            param_vals(i,:) = prior_funcs.trans_finv(param_vals_trans(i,:),args);
            parameters = param_vals(i,:);

            % simulate model and calculate discrepancy
            [param_disc(i), param_sims(i)] = args.simulate_and_calc_discrepancy(parameters,param_sims(i),[T_sim_old args.T_sim],args);
        end

        % transform the parameters
        for i=1:args.n_particles
            param_vals_trans(i,:) = prior_funcs.trans_f(param_vals(i,:),args);            %part vals is transformed
        end
    end

    % sort the particles
    [param_disc,ix] = sort(param_disc);
    param_vals_trans = param_vals_trans(ix,:);
    param_sims = param_sims(ix);

    % if most of the particles are under the final target then don't drop
    num_drop = min(sum((param_disc > dist_final)),floor(a*args.n_particles));
    num_keep = args.n_particles-num_drop;

    %calculate distances
    dist_max = param_disc(args.n_particles);
    dist_next = param_disc(num_keep);

    % check to see if we will reach desired tolerance at next iteration
    if (dist_next < dist_final)
        dist_next = dist_final;
    end
    fprintf('Simulation time=%f, max discrepancy=%f, next discrepancy=%f, number to drop=%d\n',args.T_sim,dist_max,dist_next,num_drop);

    % if acceptance is too low, then bail
    if (p_acc < p_acc_min)
        fprintf('Getting out as MCMC acceptance rate is below acceptable threshold\n');
        break;
    end

end %end SMC steps

%do the inverse transform
for i=1:args.n_particles
    param_vals(i,:) = prior_funcs.trans_finv(param_vals_trans(i,:),args);
end



end



%% MCMC function
function [params, disc,i_acc, sims] = MCMC(mcmc_steps,params,cov_matrix,args,prior_funcs,disc,i_acc,dist_next,sims)
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

    %initialise simulation storage
    prop_sims.x = 0; %saves simulation times
    prop_sims.y = prop(1:args.n_species)'; %saves simulation populations

    % simulate model and calculate discrepancy
    [dist_prop, prop_sims] = args.simulate_and_calc_discrepancy(prop,prop_sims,[0 args.T_sim],args);

    %Accept a particle if it is within the target distance.
    if (dist_prop <= dist_next)
        % then the metropolis-hastings ratio is satisfied
        params = part_vals_prop;
        disc = dist_prop;
        i_acc = i_acc + 1;
        sims = prop_sims;
    end
end

end
