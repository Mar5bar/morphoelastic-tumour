addpath(genpath('./helpers'))

%% Specify model parameters and setup. All quantities are in SI units.
params = struct();

% Model identifier, one of 'A', 'B', 'C', 'D', or 'E', corresponding to the
% accompanying publication.
params.model = 'A';

% Time and discretisation.
params.tFinal = 30; % Final time.
params.nR = 2000; % Number of points in the spatial discretisation.
params.nT = 3000; % Number of points in the temporal discretisation.

% Material parameters.
params.mu = 1; % Shear modulus of the tumour material.
params.kappa = 0; % Spring constant in the radial stress boundary condition.

% Growth rate parameters.
params.k = 1; % The basic growth rate constant.
params.sigmaHat = -1; % Threshold below which growth is arrested due to compressive stress, if included in the model.

% Nutrient parameters.
params.cInf = 1; % Concentration of nutrient at the boundary.
params.cHat = 0.8; % Nutrient threshold for necrosis.
params.lambda = 0.9; % Rate of nutrient consumption by tissue.
params.D = 1; % Diffusion coefficient of nutrient inside the tumour.
params.L = sqrt(params.D * params.cInf/params.lambda); % The diffusive lengthscale.
params.bHat = sqrt(6)*params.L; % The threshold radius of the tumour for full/partial perfusion.

% Initial configuration.
params.B = params.L; % Initial (and unstressed) radius of the spheroid.

% Timescale
params.T = 1 / (params.k * params.cInf);

% Constant thresholds for Taylor expanding integrals.
params.radialStressIntegrandThreshold = 0.05; % Radial threshold for using Taylor expansion of radial stress integrand.
params.elasticStretchIntegrandThreshold = 0.05; % Radial threshold for using Taylor expansion of of elastic stretch integrand.

%% Run the simulation.
output = runSim(params);

%% Plotting.
% plot_spheroid(output 1);
% plot_spheroid(output, params.nT);
plot_evolution(output);

%% Saving.

% Generate reduced output.
reducedOutput = struct();
reducedOutput.radii = output.rs(:,end);
reducedOutput.rFinal = output.rs(end,:);
reducedOutput.growthRatesFinal = output.growthRates(end,:);
reducedOutput.growthStretchesFinal = output.growthStretches(end,:);
reducedOutput.nutrientsFinal = output.nutrients(end,:);
reducedOutput.radialStressCentre = output.radialStresses(:,1);
reducedOutput.hoopStressCentre = output.hoopStresses(:,1);
reducedOutput.radialStressFinal = output.radialStresses(end,:);
reducedOutput.hoopStressFinal = output.hoopStresses(end,:);

% save('output.mat','output')
save('reducedOutput.mat','reducedOutput')
