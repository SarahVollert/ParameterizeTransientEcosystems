function [summ] = SaCD_feasible_equilibria(parameters,args)
%SaCD = simulate and calculate discrepancy

%% Extract information
n_species = args.n_species;

% Growth rate parameters are the first n_species parameters
r = parameters(1:n_species)';

% Reconstruct interaction matrix
A_nonzero_values = parameters(n_species+1:end);
A_values = zeros(n_species^2,1);
counter = 0;
for i=1:length(args.skip_parameters)
    if args.skip_parameters(i) ==0 %if parameter is not skipped
        counter = counter+1;
        A_values(i) = A_nonzero_values(counter);
    end
end
A = reshape(A_values,[n_species,n_species]);

%% Simulate ecosystem
% Calcuate equilibrium abundances
equilibrium_points = A\(-r);


%% Discrepancy function

%summarise in terms of feasibility and stability
measured_infeasibility = abs(sum(min(0,equilibrium_points)));
summ = measured_infeasibility;

end