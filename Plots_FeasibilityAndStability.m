function Plots_FeasibilityAndStability(samples_IC,args)

%extract info
n_samples = size(samples_IC,2);
metric_discrepancy = zeros(args.n_particles,6);  
metrics = zeros(n_samples,6); % six metrics for each sample

%calculate whether the model meets the expected behaviour
for s=1:n_samples %for each sample
    parfor p = 1:args.n_particles % for each model
        parameters = samples_IC{s}(p,:);
        parameters_noIC = parameters(args.n_species+1:end);
        metric_discrepancy(p,:) = ...
            [SaCD_feasible_equilibria(parameters_noIC,args)...              %calculate whether equilibria is feasible
            SaCD_stable_equilibria(parameters_noIC,args)...                 %calculate whether equilibria is stable
            SaCD_feasible_stable_equilibria(parameters_noIC,args)...        %calculate whether equilibria is feasible and stable
            SaCD_bounded_stable_equilibria(parameters_noIC,args)...         %calculate whether equilibria is bounded and stable
            SaCD_bounded_trajectories_SMCEEM(parameters,args)...            %calculate whether trajectory is bounded
            SaCD_bounded_slow_changing_trajectories_SMCEEM(parameters,args)]; %calculate whether trajectory is bounded and slow changing
    end
    fprintf('Finished sample %d/%d\n',s,n_samples) %print progress

    metrics(s,:) = sum(metric_discrepancy==0)./args.n_particles; %save as the percentage of particles that match the metrics
end

% Plot histogram
figure; ax = axes;
names = categorical({'(1) Feasible & stable equilibriums', ...
    '(2) Bounded & stable equilibriums', ...
    '(3) Bounded trajectories', ...
    '(4) Bounded & slow-changing trajectories'}); %axis labels for histogram 

%plot as percentages
bar(names,metrics*100)
ytickformat(ax, 'percentage');

%make colours match previous plots
CO = [linspace(0.2, 0.8,2)'.*ones(2,3); args.colors];
colororder(CO)

%label
xlabel("Constraints used to generate ensemble",'Interpreter','latex','FontSize',13)
ylabel("\begin{tabular}{c}Percentage of models in ensemble \\ with behaviours that may be expected \end{tabular}",'Interpreter','latex','FontSize',13)
legend("Feasible equilibrium", "Stable equilibrium","Feasible and stable equilibrium","Bounded and stable equilibrium","Bounded trajectory","Bounded and slow-changing trajectory",'Interpreter','latex','FontSize',10,'Location','EastOutside')
s = {'\begin{tabular}{c} (1) Positive \\ \& stable \\ equilibriums \end{tabular}',...
    '\begin{tabular}{c} (2) Bounded \\ \& stable \\ equilibriums \end{tabular}',...
    '\begin{tabular}{c} (3) Bounded \\ trajectories \end{tabular}', ...
    '\begin{tabular}{c} (4) Bounded \& \\ slow-changing \\ trajectories \end{tabular}'};
set(gca,'XTickLabel',s,'TickLabelInterpreter','latex')


end