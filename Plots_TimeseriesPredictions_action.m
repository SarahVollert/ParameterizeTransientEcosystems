function Plots_TimeseriesPredictions_action(samples_IC, args,population_multiplier,population_perturbed, T_action)

%% Calculate historic simulation
%initialise
CI = 0.95;                                                                  %95% credible interval
n_times = 100;                                                              %number of time points
times = linspace(0,args.end_time,n_times);                                  %time array
n_samples = size(samples_IC,2);                                             %ensembles to plot
EEM_abundances = zeros(args.n_species,n_times,args.n_particles,n_samples);  %for saving calculated trajectories

% set up to extract prediction medians and credible intervals
median_vals = zeros(args.n_species,n_times,n_samples);
CI_lower= zeros(args.n_species,n_times,n_samples);
CI_upper= zeros(args.n_species,n_times,n_samples);

%calculate trajectories
parfor particle = 1:args.n_particles

    for s = 1:n_samples
    %get solutions for original EEM particles
    initial_conds = samples_IC{1,s}(particle,1:args.n_species);
    parameter_vals = samples_IC{1,s}(particle,args.n_species+1:end);
    EEM_abundances(:,:,particle,s) = calculate_trajectory(parameter_vals, initial_conds, args, times);
    end

end

%find the median and credible intervals
parfor s=1:n_samples
    abund = EEM_abundances(:,:,:,s);

    %remove trajectories that couldn't be simulated
    abund = abund(:,:,not(isnan(abund(1,1,:))));

    %extract median and credible intervals
    median_vals(:,:,s) = quantile(abund,0.5,3);
    CI_upper(:,:,s) = quantile(abund,(1-CI)/2+CI,3);
    CI_lower(:,:,s) = quantile(abund,(1-CI)/2,3);
end

%% Calculate action response
%initialise
EEM_abundances_action = zeros(args.n_species,n_times,args.n_particles,n_samples);
median_vals_action = zeros(args.n_species,n_times,n_samples);
CI_lower_action= zeros(args.n_species,n_times,n_samples);
CI_upper_action= zeros(args.n_species,n_times,n_samples);
times_action = linspace(args.end_time,args.end_time+T_action,n_times);

%calculate trajectories
parfor particle = 1:args.n_particles

    for s = 1:n_samples
    %get solutions for original EEM particles
    initial_conds = EEM_abundances(:,end,particle,s)';
    initial_conds(population_perturbed) = initial_conds(population_perturbed)*population_multiplier;
    parameter_vals = samples_IC{1,s}(particle,args.n_species+1:end);
    EEM_abundances_action(:,:,particle,s) = calculate_trajectory(parameter_vals, initial_conds, args, times_action);
    end

end

%find the median and credible intervals
for s=1:n_samples
    abund_action = EEM_abundances_action(:,:,:,s);

    %remove trajectories that couldn't be simulated
    abund_action = abund_action(:,:,not(isnan(abund_action(1,1,:))));

    %extract info
    median_vals_action(:,:,s) = quantile(abund_action,0.5,3);
    CI_upper_action(:,:,s) = quantile(abund_action,(1-CI)/2+CI,3);
    CI_lower_action(:,:,s) = quantile(abund_action,(1-CI)/2,3);
end




%% Plotting one species and action

for i=1:args.n_species %repsonse species

%Figure
figure;
hold on;
tiles = tiledlayout(ceil(n_samples/2),2);
xlabel(tiles,'Time (years)','interpreter','latex')
tiles.Padding = "tight";
tiles.TileSpacing = "tight";

%find y axis limits
y_max = max([max(max(CI_upper_action(i,:,:))) max(max(CI_upper(i,:,:))) args.population_upper_bound(i)])*1.1;

    %plot results
    for s=1:n_samples
        nexttile; hold on;

        %historic
        jbfill_withLinestyle(squeeze(times),squeeze(CI_upper(i,:,s)),squeeze(CI_lower(i,:,s)),args.colors(s,:),args.colors(s,:),1,0.1,'-');
        plot(squeeze(times), squeeze(median_vals(i,:,s)),'color',args.colors(s,:),'LineWidth',1.5);

        % action
        jbfill_withLinestyle(squeeze(times_action),squeeze(CI_upper_action(i,:,s)),squeeze(CI_lower_action(i,:,s)),args.colors(s,:),args.colors(s,:),1,0.1,'--');
        plot(squeeze(times_action), squeeze(median_vals_action(i,:,s)),'color',args.colors(s,:),'LineWidth',1.5,'LineStyle','--');

        %plot bounds
        yline(args.population_lower_bound(i),'Color',[0.8 0.8 0.8],'LineWidth',1)
        yline(args.population_upper_bound(i),'Color',[0.8 0.8 0.8],'LineWidth',1)

        %plot action time
        xline(args.end_time,'--k','LineWidth',1)

        %style
        xlim([0 args.end_time+T_action])
        ylim([0 y_max])

        xlabel(args.sample_names(s),'Interpreter','latex')

%         names = [{"Hindcast ("+num2str(CI*100)+"\% credible interval)"} {'Hindcast (median)'}...
%             {"Predictions under action ("+num2str(CI*100)+"\% credible interval)"} {'Predictions under action (median)'}];
%         legend(names,'location', 'best','interpreter','latex','FontSize',10)
    end
ylabel(tiles,args.species_list(i)+" population",'interpreter','latex')

end

%% Plotting all species historically
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
        %historic
        jbfill(squeeze(times),squeeze(CI_upper(i,:,s)),squeeze(CI_lower(i,:,s)),args.colors(s,:),args.colors(s,:),1,0.1);
        plot(squeeze(times), squeeze(median_vals(i,:,s)),'color',args.colors(s,:),'LineWidth',1.5);
    end

    %plot bounds
    yline(args.population_lower_bound(i),'--k','LineWidth',1)
    yline(args.population_upper_bound(i),'--k','LineWidth',1)

    %style
    xlim([0 args.end_time])
    buffer = (args.population_upper_bound(i) - args.population_lower_bound(i))*0.15;
    ylim([args.population_lower_bound(i)-buffer args.population_upper_bound(i)+buffer])
    
    %label
    xlabel(args.species_list(i),'interpreter','latex')

    %legend
    if i==2
        names = [];
        for n=1:n_samples
            names = [names {char(args.sample_names_simple(n))}...
                {char(args.sample_names_simple(n))}];
        end
        names = [names {'Population bounds'}];

        legend(names,'location', 'northoutside','interpreter','latex','FontSize',10)
    end
end

end


