addpath(genpath('../../helpers'))
addpath(genpath('.'))

%% Run the various fitting scripts in order to produce full output.
disp('Fitting Model A to 0% gel data...')
mainModelA;
disp('Fitting Model E to 0% gel data...')
mainModelEFitTo0;
disp('Fitting Model E to 0%, 0.5% and 1% gel data...')
mainModelEPooled;

%% Package output.
disp('Packaging output...')
genModelADataForFig;
genModelEPooledDataForFig;
packageModelAParams;
packageModelEParams;
disp('Done!')