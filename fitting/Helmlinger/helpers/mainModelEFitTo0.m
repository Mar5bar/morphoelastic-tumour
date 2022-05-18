%% Initialise environment and params.
init;

% Load the optimised parameters from the Model A
% optimisation, for use as a starting guess.
outputModelA = load('outputs/outputModelA.mat');
params = outputModelA.optimParams;

% Select the correct model.
params.model = 'E';

%% Load data to fit to.
load('data/0PercentFalseInit.mat');
% Data is loaded in as radii (um).

% In this case, kappa = 0 as there is no external medium.
params.kappa = 0;

% Fit model E.
f = @(vecParams, tData) fitFuncModelETo0(vecParams, tData, params);
% Initial guess.
initGuess = [0.327      0.20347       0.3277   2.4225e-05    -0.012808];
lowerBounds = [0,0,0,0,-Inf];
upperBounds = [Inf,Inf,Inf,Inf,0];
typicalX = initGuess;

options = optimoptions('lsqcurvefit','Display','iter','TypicalX',typicalX,'UseParallel',true);
[x,resnorm,~,exitflag,output] = lsqcurvefit(f, initGuess, data.ts, data.radii, lowerBounds, upperBounds, options);

%% Plotting.
[vals, optimParams] = fitFuncModelETo0(x, data.ts, params);
output = runSim(optimParams, false);
reducedOutput = struct();
reducedOutput.ts = output.ts;
reducedOutput.rs = output.rs(:,end)*1e6;
clear output
figure
hold on
plot(data.ts,data.radii);
plot(reducedOutput.ts,reducedOutput.rs)
legend({'Data','Fit'})

%% Saving.
save('outputs/outputModelEFitTo0.mat')
