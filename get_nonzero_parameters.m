function [skip_parameters, interaction_terms_nonzero] = get_nonzero_parameters(A_sign)
% This function returns an array indicating whether a parameter in the interaction
% matrix is skipped (yes=1) or not (no=0), as well as the interaction sign
% where a term is not skipped (interaction_terms_nonzero). 

%determine which interaction terms are 0
[n_species,~] = size(A_sign);
skip_parameters = zeros(n_species^2,1);                                     %an array which indicates if an interaction term is skipped

%get all of the interaction terms as an array
interaction_terms = reshape(A_sign,[1,n_species^2]);

%extract the non-zero terms
interaction_terms_nonzero = [];
for i=1:length(interaction_terms)
    if interaction_terms(i) ==0
        skip_parameters(i)= true;
    else
        interaction_terms_nonzero = [interaction_terms_nonzero interaction_terms(i)];
    end
end

end