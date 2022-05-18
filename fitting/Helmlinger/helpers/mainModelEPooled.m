%% Initialise environment and params.
init;

% Load the optimised parameters from the fit with no external medium.
outputModelEFitTo0 = load('outputs/outputModelEFitTo0.mat');
params = outputModelEFitTo0.optimParams;
clear outputModelEFitTo0

% Select the correct model.
params.model = 'E';

%% Load data to fit to.
tData = {};
rData = {};
filenames = {'0PercentFalseInit','05Percent','1Percent'};
for i = 1 : length(filenames)
    data = load(['data/',filenames{i},'.mat']);
    tData{i} = data.data.ts;
    rData{i} = data.data.radii;
end
% Data is loaded in as radii (um).

% Fit model E.

f = @(vecParams) fitFuncModelEPooled(vecParams, tData, rData, params, 0);
% Initial guess.
% initGuess = [0.327 0.20347 0.3277 -0.012808 0.001968502946534 0.003332016445928 2.4225e-05, 2.4225e-05, 2.4225e-05];
initGuess = [0.32194      0.42638      0.26856    -0.021059    0.0021779     0.016758 4.339e-05   2.2478e-05   1.6383e-05];
lowerBounds = [0,0,0,-Inf,0,0,0,0,0];
upperBounds = [Inf,Inf,Inf,0,Inf,Inf,Inf,Inf,Inf];
typicalX = initGuess;

options = optimoptions('lsqnonlin','Display','iter','TypicalX',typicalX,'UseParallel',true);
[x,resnorm,~,exitflag,output] = lsqnonlin(f, initGuess, lowerBounds, upperBounds, options);

%% Plotting.
maxT = -Inf;
for i = 1 : length(filenames)
    maxT = max(max(tData{i}),maxT);
end
ts = linspace(0,maxT,1e3);
tDataHighRes = cell(length(filenames),1);
[tDataHighRes{:}] = deal(ts(:));
[~, rsCell] = fitFuncModelEPooled(x, tDataHighRes, cellfun(@(x)0*x,tDataHighRes,'UniformOutput',false), params);
figure
set(gcf,'Position',[90   462   945   515])
hold on
colors = {'black','blue','red'};
for i = 1 : length(filenames);
    plot(tData{i},rData{i},'o','MarkerSize',12,'Color',colors{i},'LineWidth',1,'DisplayName',['Data (',filenames{i},')']);
    plot(ts,rsCell{i},'Color',colors{i},'LineWidth',1,'DisplayName',['Fit (',filenames{i},')'])
end
legend('Location','best')
xlabel('Time (days)')
ylabel('Radius (um)')
exportgraphics(gcf,'pooledFit.png')

%% Saving.
save('outputs/outputModelEPooled.mat')