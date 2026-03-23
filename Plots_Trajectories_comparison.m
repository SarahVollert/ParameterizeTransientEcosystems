function Plots_Trajectories_comparison(samples_IC, args)
%This function plots example trajectories for each of the different sets of
%constraints in a grid.

%initialise
n_intervals = 100;                                                          %number of time intervals for calculating trajectories
n_particles_to_plot = min(size(samples_IC{1},1),30);                        %number of trajectories to plot
n_samples = size(samples_IC,2);                                             %number of constraint sets being plotted
times = linspace(0,args.end_time,n_intervals);                              %time array based on the number of intervals
EEM_abundances = zeros(args.n_species,n_intervals,args.n_particles,n_samples); %set up to save abundance info

%randomly choose parameter sets from ensemble to plot
particle_nos = randi([1,size(samples_IC{1},1)],[n_particles_to_plot,n_samples]);

%calculate trajectories
for part = 1:n_particles_to_plot %for each parameter set
    for s = 1:n_samples
        particle = particle_nos(part,s);
        %get population trajectories for parameter sets of each constraint set
        initial_conds = samples_IC{1,s}(particle,1:args.n_species);
        parameter_vals = samples_IC{1,s}(particle,args.n_species+1:end);
        EEM_abundances(:,:,part,s) = calculate_trajectory(parameter_vals, initial_conds, args, times);
    end

end

% extract population limits for setting axis limits
for sp = 1:args.n_species
    max_pops(sp) = max(max(max(EEM_abundances(sp,:,:,:))));
    max_pops(sp) = args.population_upper_bound(sp)*1.2;
    min_pops(sp) = args.population_lower_bound(sp) - (max_pops(sp) - args.population_upper_bound(sp));
end

%Figure
figure;
hold on;
tiles = tiledlayout(args.n_species,n_samples);
tiles.TileSpacing = 'compact';
tiles.Padding = 'compact';
xlabel(tiles,'Time (years)','interpreter','latex','FontSize',10)
ylabel(tiles,'Biomass/area (arbitrary units)','interpreter','latex','FontSize',10)

%for each sample
for sp=1:args.n_species

    for samp = 1:n_samples

        %for each species
        nexttile; hold on;

        %plot bounds
        xline(1,'LineStyle',':')
        area([0 1], [args.population_upper_bound(sp) args.population_upper_bound(sp)],...
            'BaseValue',args.population_lower_bound(sp), 'FaceColor','k','FaceAlpha',0.05);
        yline(args.population_lower_bound(sp),Color=[0.9 0.9 0.9],LineWidth=1.5)
        yline(args.population_upper_bound(sp),Color=[0.9 0.9 0.9],LineWidth=1.5)

        %plot trajectories
        for part = 1:n_particles_to_plot
            plot(times,EEM_abundances(sp,:,part,samp),Color=args.colors(samp,:),LineWidth=1)
        end

        %style
        xlim([0 times(end)])
        ylim([min_pops(sp) max_pops(sp)])

        %labels
        if samp == 1
            ylabel(args.species_list(sp),'interpreter','latex','rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
        end
        if sp==1
            title(args.sample_names(samp),'interpreter','latex');
        end

    end
end


end