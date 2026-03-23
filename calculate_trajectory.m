function abundances = calculate_trajectory(parameters, initial_abundances, args, times)
%This function returns the abundances for each of the species in the
%system defined by parameters and args.skip_parameters. The initial populations
%are specified by initial_abundances and the time period specified by
%priod, times is a vector of the times simulated.

%% Set up ecosystem model
% Extract information
r = parameters(1:args.n_species)'; %growth rates

%rearrange A matrix
A_nonzero_values = parameters(args.n_species+1:end);
A_values = zeros(args.n_species^2,1);
counter = 0;
for i=1:length(args.skip_parameters)
    if args.skip_parameters(i) ==0 %if parameter is not skipped
        counter = counter+1;
        A_values(i) = A_nonzero_values(counter);
    end
end
A = reshape(A_values,[args.n_species,args.n_species]);

%% Solve system 
%produce timeseries estimates
try
    solution = args.solver(@(t,n)lotkaVolterra(t,n,A,r),times,initial_abundances);
    error = false;
catch
    error = true;
end

% if the solution worked fine, get solution for specified time array (needed to compare samples at the same time point)
if error==false && solution.x(end) == times(end) %run finished
    abundances = deval(solution,times); %returns a matrix of [species x times]
else
    abundances = -ones(args.n_species,length(times))*Inf; %failed abundances
end

end