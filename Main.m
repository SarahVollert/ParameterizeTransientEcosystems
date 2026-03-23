function Main
% The following code replicates the results of Vollert, Drovandi and Adams
% (2024). For two case studies, a three-species predator-prey 
% Lotka-Volterra ecosystem model and a semiaria Australia Lotka-Volterra
% ecosystem model, this code will produce an ensemble of ecosystem models
% for each of four sets of constraints described in the manuscript. For
% trajectory-based constraints, there is an option to use SMC-EEM or the
% temporally adaptive version deveoped to speed up ensemble generation. 

% Each of the figures in the manuscript that analyse these ensembles can be
% replicated using this code. 

%% Specify case study
% Select the case study you wish to analyse from the below two options.
% 1 = 3 species case study
% 2 = semiarid case study

CASE_STUDY = 1; 

if CASE_STUDY ==1
    %set up for the 3 species case study
    excel_file_name = 'interactions_3species.xlsx';                         % This file contains the interaction network. 
    args.population_upper_bound = [15; 25; 30];                             % User specified population bounds for each species
    args.population_lower_bound = [5; 0; 0];                                % User specified population bounds for each species
    args.T_max = 1;                                                         % the observation period for trajectory-based constraints
    args.R_max = [1; 2; 6];                                                 % User specified limits of population increase rate (per one time step)
    args.R_min = [-1; -2; -6];                                              % User specified limits of population decrease rate (per one time step)
    args.n_particles = 1000;                                               % Ensemble size
    args.species_list_simple = [{'A'} {'B'} {'C'}];                         % Simple species labels for figures

elseif CASE_STUDY ==2 
    %set up for the 3 species case study
    excel_file_name = 'interactions_SemiaridAustralia.xlsx';                % This file contains the interaction network. 
    args.population_upper_bound = [100; 500; 600; 1000; 1000; 1000; 10; 1000]; % User specified population bounds for each species
    args.population_lower_bound = [10; 10; 10; 50; 100; 100; 0; 0];         % User specified population bounds for each species
    args.T_max = 1;                                                         % the observation period for trajectory-based constraints
    args.R_max = 12*ones(size(args.population_upper_bound));                % User specified limits of population increase rate (per one time step)
    args.R_min = -36*ones(size(args.population_upper_bound));               % User specified limits of population decrease rate (per one time step)
    args.n_particles = 50000;                                               % Ensemble size
    args.species_list_simple = [{'D'} {'M'} {'LH'} {'SV'} {'G'} {'I'} {'F'} {'S'}]; % Simple species labels for figures 
end

%% Specify ecosystem information
args.upper_growth_rate = 5;                                                 % maximum value of r for each species
args.lower_growth_rate = -5;                                                % minimum value of r for each species
args.upper_interaction_strength = 1;                                        % maximum absolute value of alpha for each interaction
args.solver = @ode45;                                                       % ODE solver can be specified differently

%% Specify search information
prior_funcs = define_prior_uniform;                                         % defines the prior distribution as uniform

%% Prepare to run search
%extract information for this ecosystem
args = prepare_args_for_sampling(args, excel_file_name);                    % sets up the ecosystem model as specified by interaction network excel file                 

%% Specify plotting details
[~,args.species_list] = xlsread(excel_file_name,'A:A');                     %define the names of the species for plotting
args.param_names = define_parameter_names(args.skip_parameters,args.species_list_simple); %define parameter names for parameter plots
args.max_per_figure = args.num_params;                                      %the maximum number of parameters to include in a single plot
args.text_not_label = false;                                                %if true, add text to figure, if false species indicated by label 
args.end_time = 1;                                                          %the period for timeseries plotting
args.colors = [0.8510    0.3255    0.0980;                                  %colour definitions for each discrepancy option
                0.9608    0.7529    0.2706; 
                0.4667    0.6745    0.1882; 
                0.0745    0.6235    1.0000];
args.ICs = rand(args.n_particles,args.n_species).*(args.population_upper_bound-args.population_lower_bound)' +args.population_lower_bound'; % random initial conditions for plotting equilibrium constraints

%% GENERATING AN ENSEMBLE OF PARAMETER SETS FOR EACH CONSTAINT SET
% 1. feasible and stable equilibria
args.simulate_and_calc_discrepancy = @SaCD_feasible_stable_equilibria;      %select the appropriate constaint set for simulation and discrepancy calculation
tic;[EEM_sample_1, prior_sample] = SMCEEM_method(args,prior_funcs);toc %sample using SMC-EEM and time
EEM_sample_1_IC = [args.ICs EEM_sample_1];                    %save a version with random initial conditions attached

% 2. bounded and stable equilibria
args.simulate_and_calc_discrepancy = @SaCD_bounded_stable_equilibria; %select the appropriate constaint set for simulation and discrepancy calculation
tic; [EEM_sample_2, ~] = SMCEEM_method(args,prior_funcs);toc       %sample using SMC-EEM and time
EEM_sample_2_IC = [args.ICs EEM_sample_2];                %save a version with random initial conditions attached

% 3. bounded trajectories
%Update prior bounds to include the initial conditions
args.lower = [args.population_lower_bound' args.lower];    %lower lim
args.upper = [args.population_upper_bound' args.upper];    %upper lim
args.num_params = size(args.upper,2);
    
    %SMC-EEM
    args.simulate_and_calc_discrepancy = @SaCD_bounded_trajectories_SMCEEM;  %select the appropriate constaint set for simulation and discrepancy calculation
    tic;[EEM_sample_3_IC_SMCEEM, prior_sample_IC] = SMCEEM_method(args,prior_funcs);toc %sample using SMC-EEM and time
    EEM_sample_3_SMCEEM = EEM_sample_3_IC_SMCEEM(:,args.n_species+1:end);  %save a version without initial conditions
    
    % SMC w/ adaptive time stepping 
    args.simulate_and_calc_discrepancy = @SaCD_bounded_trajectories_TASMCEEM; %select the appropriate constaint set for simulation and discrepancy calculation
    tic;[EEM_sample_3_IC_TASMCEEM, ~] = SMCEEM_method_temporallyAdaptive(args,prior_funcs); toc; %sample using temporally adaptive SMC-EEM and time
    EEM_sample_3_TASMCEEM = EEM_sample_3_IC_TASMCEEM(:,args.n_species+1:end);  %save a version without initial conditions

% 4. bounded and slow-changing trajectories
    %SMC-EEM
    args.simulate_and_calc_discrepancy = @SaCD_bounded_slow_changing_trajectories_SMCEEM; %select the appropriate constaint set for simulation and discrepancy calculation
    tic;[EEM_sample_4_IC_SMCEEM, ~] = SMCEEM_method(args,prior_funcs); toc  %sample using SMC-EEM and time
    EEM_sample_4_SMCEEM = EEM_sample_4_IC_SMCEEM(:,args.n_species+1:end);  %save a version without initial conditions

    % SMC w/ adaptive time stepping 
    args.simulate_and_calc_discrepancy = @SaCD_bounded_slow_changing_trajectories_TASMCEEM; %select the appropriate constaint set for simulation and discrepancy calculation
    tic;[EEM_sample_4_IC_TASMCEEM, ~] = SMCEEM_method_temporallyAdaptive(args,prior_funcs);toc %sample using temporally adaptive SMC-EEM and time
    EEM_sample_4_TASMCEEM = EEM_sample_4_IC_TASMCEEM(:,args.n_species+1:end);  %save a version without initial conditions

%% LOAD IN SEMIARID CASE STUDY RESULTS FROM MANUSCRIPT
if CASE_STUDY == 1
      load('RESULTS_3Species.mat');
else  load('RESULTS_SemiaridAustralia.mat');
end

% Set up constraint comparison samples for plotting
args.sample_names = [{'\begin{tabular}{c} (1) Positive \\ \& stable \\ equilibriums \end{tabular}'} ...
    {'\begin{tabular}{c} (2) Bounded \\ \& stable \\ equilibriums \end{tabular}'} ...
    {'\begin{tabular}{c} (3) Bounded \\ trajectories \end{tabular}'} ...
    {'\begin{tabular}{c} (4) Bounded \& \\ slow-changing \\ trajectories \end{tabular}'}];
args.sample_names_simple = [{'(1) Positive \& stable equilibriums'} ...
    {'(2) Bounded \& stable equilibriums'} ...
    {'(3) Bounded trajectories'} ...
    {'(4) Bounded \& slow-changing trajectories'}];
samples = [{EEM_sample_1_IC(:,args.n_species+1:end)}...
    {EEM_sample_2_IC(:,args.n_species+1:end)} ...
    {EEM_sample_3_IC_SMCEEM(:,args.n_species+1:end)} ...
    {EEM_sample_4_IC_SMCEEM(:,args.n_species+1:end)}];
samples_IC = [{EEM_sample_1_IC} {EEM_sample_2_IC} {EEM_sample_3_IC_SMCEEM} {EEM_sample_4_IC_SMCEEM}];



%% Plot comparison of constraints
% Individual parameter trajectories
Plots_Trajectories_comparison(samples_IC, args); 

% Plot parameter distributions
args.IC_plotting = true;
Plots_ParameterDistributions_comparison(prior_sample_IC, samples_IC, args);

% Model hindcast + predict response to a conservation action
    %action inputs
    population_multiplier = 0.5;                                            %percentage of population remaining at pulse perturbation
    population_perturbed = 1;                                               %the species number to be perturbed
    args.end_time = args.T_max;                                             %historical period to model before action
    T_action = 1/2;                                                         % time to model following action
Plots_TimeseriesPredictions_action(samples_IC, args,population_multiplier,population_perturbed, T_action)

% Plot percentage feasible and stable in ensembles
Plots_FeasibilityAndStability(samples_IC,args)

%% Plot comparison of algorithms
% Set up algorithm comparisons
args.colors = [0 0 0; 0 0 1];
args.sample_names_simple = [{'(1) SMC-EEM'} ...
    {'(2) Temporally adaptive SMC-EEM'} ];
samples_IC_algorithms = [{EEM_sample_3_IC_SMCEEM} {EEM_sample_3_IC_TASMCEEM}];

%Plots
Plots_TimeseriesPredictions_comparison_IC(samples_IC_algorithms,args); 
Plots_ParameterDistributions_comparison(prior_sample_IC, samples_IC_algorithms, args); 

end