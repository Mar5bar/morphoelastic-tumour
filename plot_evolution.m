function plot_evolution(ts, rs, growthStretches, growthRates, stresses, nutrients, necroticRadii, params)
%% Plot the evolution of various scalar quantities over time.

thetas = linspace(0,pi/2,1e2)';
x = @(r) r.*cos(thetas);
y = @(r) r.*sin(thetas);


figure
t = tiledlayout('flow');

% Outer radius and necrotic radius.
nexttile()
hold on
plot(ts,rs(:,end),'Color','black','LineWidth',1)
plot(ts, necroticRadii,'Color',0.7*[1,1,1],'LineWidth',1);
legendEntries = {'Outer radius','Necrotic radius'};
if any(nutrients(:) == 0)
    nutrientsTemp = nutrients; nutrientsTemp(nutrients > 0) = 1;
    [~,zeroNutrientThresholdInds] = max(nutrientsTemp,[],2,'linear');
    plot(ts, rs(zeroNutrientThresholdInds),'Color',0.4*[1,1,1],'LineWidth',1);
    legendEntries{end+1} = 'Nutrient-free radius';
end
legend(legendEntries,'Location','southeast')
box on
xlabel('$t$')
title('Spheroid radius')