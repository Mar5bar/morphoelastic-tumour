function plot_spheroid(ts, rs, growthStretches, growthRates, stresses, nutrients, params, frame)
%% Plot a quarter spheroid with various shadings.
if nargin < 8
    frame = size(rs,1);
end

thetas = linspace(0,pi/2,1e2)';
x = @(r) r.*cos(thetas);
y = @(r) r.*sin(thetas);

shaders = {growthStretches, growthRates, stresses, nutrients};
shaderStrings = {
            '$\gamma$',...
            '$\frac{1}{\gamma}\frac{\partial\gamma}{\partial t}$',...
             '$\sigma_r$',...
             '$c$'...
             };
thresholdInds = {
            find(growthStretches(frame,:) >= 1,1,'first'),...
            find(growthRates(frame,:) >= 0,1,'first'),...
            find(stresses(frame,:) >= params.stressGrowthThreshold,1,'first'),...
            find(nutrients(frame,:) >= params.necrosisThreshold,1,'first')
            };
thresholdLabels = {
            '$1$',...
            '$0$',...
            '$\hat{\sigma}_r$',...
            '$\hat{c}$'...
            };

figure
t = tiledlayout('flow');
for i = 1 : length(shaders)

    nexttile()
    hold on
    xs = x(rs(frame,:)); ys = y(rs(frame,:)); vals = repmat(shaders{i}(frame,:),length(thetas),1);
    mx = max(xs(:)); my = max(ys(:));
    % Shade.
    pcolor(xs,ys,vals)
    shading interp
    colormap(viridis)
    c = colorbar;
    set(c,'TickLabelInterpreter','latex')
    % Outer radius.
    plot(x(rs(frame,end)),y(rs(frame,end)),'Color','black','LineWidth',1)
    if thresholdInds{i}
        xs = x(rs(frame,thresholdInds{i})); ys = y(rs(frame,thresholdInds{i}));
        plot(xs,ys,'Color','black','LineWidth',1)
        xLoc = xs(end/2)-mx/20;
        yLoc = ys(end/2)-my/20;
        % Only plot a label if it would be visible.
        if xLoc>0 & yLoc>0 & thresholdInds{i}<size(rs,2)
            text(xLoc,yLoc,thresholdLabels{i},'Interpreter','latex','FontSize',20)
        end
    end
    title(shaderStrings{i})
    axis equal
    axis tight

end
title(t,['Spheroid at $t=',num2str(ts(frame)),'$'],'Interpreter','latex','FontSize',24)
set(gcf,'Position',[391   210   864   767])