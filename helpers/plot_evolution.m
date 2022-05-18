function plot_evolution(output)
%% Plot the evolution of various scalar quantities over time.

thetas = linspace(0,pi/2,1e2)';
x = @(r) r.*cos(thetas);
y = @(r) r.*sin(thetas);


figure
t = tiledlayout('flow');

% Outer radius and necrotic radius.
nexttile()
hold on
plot(output.ts / output.params.T,output.rs(:,end) / output.params.L,'Color','black','LineWidth',1)
plot(output.ts / output.params.T, output.necroticRadii / output.params.L,'Color',0.7*[1,1,1],'LineWidth',1);
legendEntries = {'Outer radius','Necrotic radius'};
if any(output.nutrients(:) == 0)
    nutrientsTemp = output.nutrients; nutrientsTemp(output.nutrients > 0) = 1;
    [~,zeroNutrientThresholdInds] = max(nutrientsTemp,[],2,'linear');
    plot(output.ts, output.rs(zeroNutrientThresholdInds) / output.params.L,'Color',0.4*[1,1,1],'LineWidth',1);
    legendEntries{end+1} = 'Nutrient-free radius';
end
legend(legendEntries,'Location','southeast')
box on
xlabel('$t/T$')
ylabel('$r/L$')
title('Spheroid radius')