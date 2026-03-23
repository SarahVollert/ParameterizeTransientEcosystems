function prior_funcs = define_prior_uniform
% This function defines all of the functions needed to specify the prior
% distribution and returns them all in the prior_funcs structure. 

%define prior
    prior_funcs.sampler = @(args) ...
        [unifrnd(args.lower,args.upper)];                       % uniform distribution
    %define transform
    prior_funcs.trans_f = @(theta,args) ...
        [uniform_transform(theta,args)];                              % uniform transform
    prior_funcs.trans_finv = @(theta_trans,args) ...
        uniform_transform_inverse(theta_trans,args);                  % uniform transform
    %define the prior pdf
    prior_funcs.pdf = @(theta_trans,args) ...                         % uniform transformed pdf
        uniform_pdf_transformed(theta_trans);

end


function theta = uniform_pdf_transformed(theta_trans)
    % This function calculates the density of the transformed uniform
    % distribution
        theta = prod(abs(exp(theta_trans)./(exp(theta_trans)+1).^2));
end


function theta_transformed = uniform_transform(theta,args)
    % This function uses a transform for uniform distributions, so that all
    % values are within the uniform bounds.
    theta_transformed = zeros(size(theta));
    for i=1:size(theta,2)
        theta_transformed(:,i) = log((theta(:,i)-args.lower(i))./(args.upper(i)-theta(:,i)));
    end
end


function theta = uniform_transform_inverse(theta_trans,args)
    % This function inverses the transform for uniform distributions.
    theta = zeros(size(theta_trans));
    for i=1:size(theta_trans,2)
        theta(:,i) = (args.lower(i)+args.upper(i)*exp(theta_trans(:,i)))./(1+exp(theta_trans(:,i)));
    end
end