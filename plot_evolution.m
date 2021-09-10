function plot_evolution(ts, rs, growthStretches, growthRates, stresses, nutrients, params)
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
nutrientsTemp = nutrients; nutrientsTemp(nutrients<params.necrosisThreshold) = 0; nutrientsTemp(nutrients>=params.necrosisThreshold) = 1;
[~,necroticThresholdInds] = max(nutrientsTemp,[],2,'linear');
plot(ts, rs(necroticThresholdInds),'Color',0.7*[1,1,1],'LineWidth',1);
legend({'Outer radius','Necrotic radius'},'Location','southeast')
box on
xlabel('$t$')
title('Spheroid radius')