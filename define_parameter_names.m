function param_names = define_parameter_names(skip_parameters,species_list)
% This function defines the parameter names based on which parameters
% remain in the model. 

%set up
n_species = sqrt(length(skip_parameters));

%growth rates 
for i=1:n_species
    param_names(i) = "$r_{"+species_list(i)+"}$";                           %set up to be plotted using latex interpreter
end

%interaction terms
skip_parameters_matrix = reshape(skip_parameters,[n_species,n_species]);
counter = 0;
for i=1:n_species
    for j=1:n_species
        if skip_parameters_matrix(i,j) ==0
            counter = counter+1;
            param_names(n_species+counter) = "$\alpha_{"+species_list(j)+","+species_list(i)+"}$"; %set up to be plotted using latex interpreter
        end
    end
end

end