%% Initialise environment and params.
init;

% Load the optimised parameters from the fit with no external medium.
outputModelEFitTo0 = load('outputs/outputModelEFitTo0.mat');
params = outputModelEFitTo0.optimParams;
clear outputModelEFitTo0

% Select the correct model.
params.model = 'E';

%% Load data to fit to.
load('data/1Percent.mat');
% Data is loaded in as radii (um).

% Fit model E.
f = @(vecParams, tData) fitFuncModelE(vecParams, tData, params);
% Initial guess.
initGuess = [0.001541921985180];
lowerBounds = [0];
upperBounds = [1];
typicalX = initGuess;

options = optimoptions('lsqcurvefit','Display','iter','TypicalX',typicalX,'UseParallel',true);
[x,resnorm,~,exitflag,output] = lsqcurvefit(f, initGuess, data.ts, data.radii, lowerBounds, upperBounds, options);

%% Plotting.
[vals, optimParams] = fitFuncModelE(x, data.ts, params);
output = runSim(optimParams, false);
figure
hold on
plot(data.ts,data.radii,'o','MarkerSize',12,'Color','black','LineWidth',1);
plot(output.ts,output.rs(:,end)*1e6,'Color','black','LineWidth',1)
legend({'Data','Fit'})

%% Saving.
save('outputs/outputModelE.mat')
