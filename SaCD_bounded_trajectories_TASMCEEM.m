function [discrepancy, previous_sims_updated] = SaCD_bounded_trajectories_TASMCEEM(parameters,previous_sims,period,args)
%SaCD = simulate and calculate discrepancy
% Bounded populations obtained numerically, rather than by calculating 
% equilibium populations. This function simulates the ecosystem features 
% and summarises them in the discrepancy function. 

%This version is suitable for the temporally adaptive SMC-EEM algorithm (TASMCEEM)

%% Extract information
n_species = args.n_species;

% Growth rate parameters are the second n_species parameters
r = parameters(n_species+1:n_species*2)';

% Reconstruct interaction matrix
A_nonzero_values = parameters(n_species*2+1:end);
A_values = zeros(n_species^2,1);
counter = 0;
for i=1:length(args.skip_parameters)
    if args.skip_parameters(i) ==0 %if parameter is not skipped
        counter = counter+1;
        A_values(i) = A_nonzero_values(counter);
    end
end
A = reshape(A_values,[n_species,n_species]);

%% Solve system
warning('off','all')
time_array = linspace(period(1), period(2), 100);
Opt = odeset('Events', @explosivePopulationEvent);

%solve system from initial conditions of starting_pops for time specified by period
try
%     solution = args.solver(@(t,n) lotkaVolterra(t,n,A,r), time_array, previous_sims.y(:,end));
    solution = args.solver(@(t,n,args) lotkaVolterra(t,n,A,r), time_array, previous_sims.y(:,end), Opt, args);
    error = false;

    %partial solutions cannot be treated as full solutions
    if any([any(isnan(solution.y)) solution.x(end)<period(end)])
        error = true;
    end 

catch
    error = true;
end

%% If there are errors
if error == true
    %update simulation
    previous_sims_updated.x = [previous_sims.x previous_sims.x(end)*1.01]; 
    previous_sims_updated.y = [previous_sims.y args.population_upper_bound*Inf];

    %set discrepancy
    discrepancy = Inf;
else
    %update simulation
    previous_sims_updated.x = [previous_sims.x solution.x(2:end)]; %go from 2 to avoid double records
    previous_sims_updated.y = [previous_sims.y solution.y(:,2:end)];

    %calculate total area above max population line
    populations_above = max(previous_sims_updated.y - args.population_upper_bound,0);
    for sp = 1:args.n_species
        area_above(sp) = trapz(previous_sims_updated.x, populations_above(sp,:));
    end

    %calculate total area below min population line
    populations_below = max(args.population_lower_bound-previous_sims_updated.y,0);
    for sp = 1:args.n_species
        area_below(sp) = trapz(previous_sims_updated.x, populations_below(sp,:));
    end

    measured_population_deficiency = sum(area_below);
    measured_population_surplus = sum(area_above);
    discrepancy = measured_population_deficiency+ measured_population_surplus;
end
end


