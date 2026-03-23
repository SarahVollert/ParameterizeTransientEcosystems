function [value, isterminal, direction] = explosivePopulationEvent(~,n,args)
% This is an event function to sit within an ODE solver which stops the
% simulation from running if explosive populations are encountered. 

% Explosive populations are determined by whether the poopulation upper
% bound is exceeded by 4 orders of magnitude for any species


% thresholds for which we care about (each value(i) describes the ith event)
value(1:args.n_species) = args.population_upper_bound'*10^4-n'; %pops too high

% stopping conditions
isterminal = ones(1,args.n_species); %1 indicates stop the solve when the ith event occues
direction = zeros(1,args.n_species); %0 = can hit threshold from either direction (1= only if function is increasing to value, -1 if decreasing to value)
end
