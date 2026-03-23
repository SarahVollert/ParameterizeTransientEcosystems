function Plots_TimeseriesPredictions_comparison_IC(samples_IC, args)
% This function plots the timeseries predictions for any samples passed.

%set credible interval
CI =0.95;

%initialise
n_intervals = 1000;
times = linspace(0,args.end_time,n_intervals);
n_samples = size(samples_IC,2);
EEM_abundances = zeros(args.n_species,n_intervals,args.n_particles,n_samples);

%calculate trajectories
for particle = 1:args.n_particles

    for s = 1:n_samples
    %get solutions for original EEM particles
    initial_conds = samples_IC{1,s}(particle,1:args.n_species);
    parameter_vals = samples_IC{1,s}(particle,args.n_species+1:end);
    EEM_abundances(:,:,particle,s) = calculate_trajectory(parameter_vals, initial_conds, args, times);
    end

end

%find the median and credible intervals
for s=1:n_samples
    abund = EEM_abundances(:,:,:,s);

    %remove trajectories that couldn't be simulated
    abund = abund(:,:,not(isnan(abund(1,1,:))));

    %extract info
    median_vals(:,:,s) = quantile(abund,0.5,3);
    CI_lower(:,:,s) = quantile(abund,(1-CI)/2+CI,3);
    CI_upper(:,:,s) = quantile(abund,(1-CI)/2,3);
end

%%
%Figure
figure;
hold on;
tiles = tiledlayout(ceil(args.n_species/3),3);
xlabel(tiles,'Time (years)','interpreter','latex')
ylabel(tiles,'Populations','interpreter','latex')

%for each species
for i=1:args.n_species
    nexttile; hold on;

    %plot results
    for s=1:n_samples
        if s==1
            jbfill_withLinestyle(squeeze(times),squeeze(CI_upper(i,:,s)),squeeze(CI_lower(i,:,s)),args.colors(s,:),args.colors(s,:),1,0.1,'-');
            plot(squeeze(times), squeeze(median_vals(i,:,s)),'color',args.colors(s,:),'LineWidth',1.5);
        else
            jbfill_withLinestyle(squeeze(times),squeeze(CI_upper(i,:,s)),squeeze(CI_lower(i,:,s)),args.colors(s,:),args.colors(s,:),1,0.1,'--');
            plot(squeeze(times), squeeze(median_vals(i,:,s)),'color',args.colors(s,:),'LineWidth',1.5,'LineStyle','--');
        end
    end

    %plot bounds
    yline(args.population_lower_bound(i),'--k','LineWidth',1)
    yline(args.population_upper_bound(i),'--k','LineWidth',1)


    %style
    xlim([0 args.end_time])
    
    %label
    xlabel(args.species_list(i),'interpreter','latex')

    %legend
    if i==2
        int1_text = " ("+num2str(CI*100)+"\% credible interval)";
        names = [];
        for n=1:n_samples
            names = [names {char(args.sample_names(n))+int1_text}...
                {char(args.sample_names(n))+" (median)"}];
        end
        names = [names {'Population bounds'}];

        legend(names,'location', 'northoutside','interpreter','latex','FontSize',10)
    end
end


end


