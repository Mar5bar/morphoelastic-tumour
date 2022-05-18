loaded = load('outputs/outputModelA.mat');
params = loaded.optimParams;
save('fittedParameters/ModelAFittedParameters-0%.mat','params')