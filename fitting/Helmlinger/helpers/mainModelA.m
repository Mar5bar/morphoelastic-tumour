%% Initialise environment and params.
init;

%% Load data to fit to.
load('data/0PercentFalseInit.mat');
% Data is loaded in as radii (um).

% Select the correct model.
params.model = 'A';

% Fit a basic Model A.
f = @(vecParams, tData) fitFuncModelA(vecParams, tData, params);
% Initial guess.
initGuess = [0.38726      0.13956      0.81783    1.833e-05];
lowerBounds = [0,0,0,0];
upperBounds = [Inf,Inf,Inf,Inf];
typicalX = initGuess;

options = optimoptions('lsqcurvefit','Display','iter','TypicalX',typicalX,'UseParallel',true);
[x,resnorm,~,exitflag,output] = lsqcurvefit(f, initGuess, data.ts, data.radii, lowerBounds, upperBounds, options);

%% Plotting.
[vals, optimParams] = fitFuncModelA(x, data.ts, params);
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
save('outputs/outputModelA.mat')