function [summ] = SaCD_bounded_stable_equilibria(parameters,args)
%SaCD = simulate and calculate discrepancy
% This function simulates the ecosystem features and summarises them in the
% discrepancy function for bounded and stable ecosystems. 
% 
% To simulate the ecosystem, the equilibrium abundances are calulated, and 
% then also used to calculate the eigenvalues of the Jacobian matrix. This 
% information is then summarised using the discrepancy function 
% (function output 'summ'). 

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

% Calculate Jacobian matrix
jacobian = zeros(n_species,n_species);
for k = 1:n_species
    jacobian(k,:) = A(k,:)*equilibrium_points(k);
    jacobian(k,k) = jacobian(k,k) + sum(A(k,:).*equilibrium_points');
end
jacobian(eye(n_species)==1) = jacobian(eye(n_species)==1) + r;

% Calculate real parts of the eigenvalues
real_eigenvalues = real(eig(jacobian));

%% Discrepancy function

%summarise in terms of feasibility and stability
measured_infeasibility = sum(max(0,equilibrium_points-args.population_upper_bound))+ sum(max(0,args.population_lower_bound-equilibrium_points));
measured_instability = sum(max(0,real_eigenvalues));

summ = measured_infeasibility+measured_instability;

end
