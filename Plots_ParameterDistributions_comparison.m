function Plots_ParameterDistributions_comparison(prior_sample, samples, args)
% This script plots the estimated marginal posterior densities for both
% all methods entered in samples


%figure layout settings
n_columns = round(sqrt(args.max_per_figure)/2)*2+1; %make an odd number, roughly square
n_rows = ceil(args.max_per_figure/n_columns);
n_samples = size(samples,2);
[~, n_params] = size(samples{1,1});

if args.IC_plotting == true
    for i=1:args.n_species
        IC_labels(i) = "$n_{"+args.species_list_simple(i)+"}(0)$";
    end
    full_param_names = [IC_labels args.param_names];
end

%for each parameter
for i=1:n_params

    %Do we need a new figure?
    if(mod(i,args.max_per_figure)==1)
        %set up new figure
        figure;
        tiles = tiledlayout(n_rows,n_columns);
        tiles.Padding = 'tight';
    end
    
    %begin plotting
    nexttile
    hold on

    %axis labels and limits
    set(gca,'ytick',[])
    if args.IC_plotting == true
        xlim([args.lower(i) args.upper(i)]);
        xlabel(full_param_names(i),'Interpreter','latex','FontSize',12)
    else
        xlim([args.lower_noIC(i) args.upper_noIC(i)]);
        xlabel(args.param_names(i),'Interpreter','latex','FontSize',12)
    end

    
    %plot the prior distribution in grey
    if args.IC_plotting == true
        [priorPlot_f,priorPlot_x] = ksdensity(prior_sample(:,i),'BoundaryCorrection','reflect','Support',[args.lower(i) args.upper(i)]);
    else
        [priorPlot_f,priorPlot_x] = ksdensity(prior_sample(:,i),'BoundaryCorrection','reflect','Support',[args.lower_noIC(i) args.upper_noIC(i)]);
    end
    priorPlot_f(1) = priorPlot_f(2);priorPlot_f(end)=priorPlot_f(end-1);
    area(priorPlot_x,priorPlot_f,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);


    %plot each of the samples 
    for j=1:n_samples
        if args.IC_plotting == true
            [EEMPlot_f, EEMPlot_x] = ksdensity(samples{1,j}(:,i),'BoundaryCorrection','reflection','Support',[args.lower(i) args.upper(i)]);
        else
            [EEMPlot_f, EEMPlot_x] = ksdensity(samples{1,j}(:,i),'BoundaryCorrection','reflection','Support',[args.lower_noIC(i) args.upper_noIC(i)]);
        end        
        EEMPlot_f(1) = EEMPlot_f(2);EEMPlot_f(end)=EEMPlot_f(end-1);
        plot(EEMPlot_x,EEMPlot_f,'Color',args.colors(j,:),'LineWidth',1)
    end
    
    %add legend
    if i==ceil(n_columns/2)
        names = [{'Prior sample'} args.sample_names_simple];
        legend(names(:),'Location','northoutside','Interpreter','latex','FontSize',11,'Orientation','vertical');
    end
end

%y axis label
ylabel(tiles,'Density (Arbitrary Units)','Interpreter','latex','FontSize',12)
xlabel(tiles,'Parameters','Interpreter','latex','FontSize',12)

end