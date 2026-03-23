function args = prepare_args_for_sampling(args, excel_file_name)
% This function sets up the modelling problem for the ecosystem network
% specified by the interaction matrix excel file. 

%Extract information from excel
[args.A_sign,~,~] = xlsread(excel_file_name);                                    %interaction signs matrix, A_sign
[args.n_species,~] = size(args.A_sign);                                          %number of species in ecosystem network
[args.skip_parameters, args.interaction_terms_nonzero] = get_nonzero_parameters(args.A_sign);                 % the parmaeters that are zero as defined by A_sign

%interaction bounds - based on interaction signs
lower_interaction_bound = args.interaction_terms_nonzero;
upper_interaction_bound = args.interaction_terms_nonzero;
for i=1:length(args.interaction_terms_nonzero)
    if args.interaction_terms_nonzero(i) == 1
        lower_interaction_bound(i)=0;
    elseif args.interaction_terms_nonzero(i) == -1
        upper_interaction_bound(i)=0;
    end
end

%prior bounds
args.lower = [ones(1,args.n_species)*args.lower_growth_rate ...
    lower_interaction_bound*args.upper_interaction_strength];    %lower lim
args.upper = [ones(1,args.n_species)*args.upper_growth_rate ...
    upper_interaction_bound*args.upper_interaction_strength];    %upper lim

% save info to pass
args.num_params = args.n_species+args.n_species^2 - sum(args.skip_parameters);       %number of parameters

% this is the lower and upper bounds for each parameter, without including
% initial conditions
args.lower_noIC = args.lower;    %lower lim
args.upper_noIC = args.upper;    %upper lim

end