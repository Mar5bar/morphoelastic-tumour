function plot_spheroid(output,frame)
%% Plot a quarter spheroid with various shadings.
if nargin < 2
    frame = size(output.rs,1);
end

thetas = linspace(0,pi/2,1e2)';
x = @(r) r.*cos(thetas);
y = @(r) r.*sin(thetas);

shaders = {output.growthStretches, output.growthRates, output.radialStresses, output.nutrients};
shaderStrings = {
            '$\gamma$',...
            '$\frac{1}{\gamma}\frac{\partial\gamma}{\partial t}$',...
             '$\sigma_r$',...
             '$c$'...
             };
thresholdInds = {
            find(output.growthStretches(frame,:) >= 1,1,'first'),...
            find(output.growthRates(frame,:) >= 0,1,'first'),...
            find(output.radialStresses(frame,:) >= output.params.sigmaHat,1,'first'),...
            find(output.nutrients(frame,:) >= output.params.cHat,1,'first')
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
    xs = x(output.rs(frame,:)); ys = y(output.rs(frame,:)); vals = repmat(shaders{i}(frame,:),length(thetas),1);
    mx = max(xs(:)); my = max(ys(:));
    % Shade.
    pcolor(xs,ys,vals)
    shading interp
    colormap(viridis)
    c = colorbar;
    set(c,'TickLabelInterpreter','latex')
    % Outer radius.
    plot(x(output.rs(frame,end)),y(output.rs(frame,end)),'Color','black','LineWidth',1)
    if thresholdInds{i}
        xs = x(output.rs(frame,thresholdInds{i})); ys = y(output.rs(frame,thresholdInds{i}));
        plot(xs,ys,'Color','black','LineWidth',1)
        xLoc = xs(end/2)-mx/20;
        yLoc = ys(end/2)-my/20;
        % Only plot a label if it would be visible.
        if xLoc>0 & yLoc>0 & thresholdInds{i}<size(output.rs,2)
            text(xLoc,yLoc,thresholdLabels{i},'Interpreter','latex','FontSize',20)
        end
    end
    title(shaderStrings{i})
    axis equal
    axis tight

end
title(t,['Spheroid at $t=',num2str(output.ts(frame)),'$'],'Interpreter','latex','FontSize',24)
set(gcf,'Position',[391   210   864   767])