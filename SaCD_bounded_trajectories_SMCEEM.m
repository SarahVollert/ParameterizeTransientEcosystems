function [summ] = SaCD_bounded_trajectories_SMCEEM(parameters,args)
%SaCD = simulate and calculate discrepancy
% Bounded populations obtained numerically, rather than by calculating 
% equilibium populations. This function simulates the ecosystem features 
% and summarises them in the discrepancy function. 

%This version is NOT suitable for the temporally adaptive SMC-EEM algorithm

%% Extract information
n_species = args.n_species;

% Initial condition parameters are the first n_species parameters
IC = parameters(1:n_species)';

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
time_array = linspace(0, args.T_max, 10000);
Opt = odeset('Events', @explosivePopulationEvent);

%solve system from initial conditions
try
%     solution = args.solver(@(t,n) lotkaVolterra(t,n,A,r), time_array, IC);
    solution = args.solver(@(t,n,args) lotkaVolterra(t,n,A,r), time_array, IC, Opt, args);
    error = false;

    %partial solutions cannot be treated as full solutions
    if any([any(isnan(solution.y)) solution.x(end)<args.T_max])
        error = true;
    end
catch
    error = true;
end


%% Measure discrepancy
if error ==false

    %calculate total area above max population line
    populations_above = max(solution.y - args.population_upper_bound,0);
    for sp = 1:args.n_species
        area_above(sp) = trapz(solution.x, populations_above(sp,:));
    end

    %calculate total area below min population line
    populations_below = max(args.population_lower_bound-solution.y,0);
    for sp = 1:args.n_species
        area_below(sp) = trapz(solution.x, populations_below(sp,:));
    end 

    %discrepancy measure
    measured_population_deficiency = sum(area_below);
    measured_population_surplus = sum(area_above);
    summ = measured_population_deficiency+ measured_population_surplus;

else
    summ = Inf; 
end


end


