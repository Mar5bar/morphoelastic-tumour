clear all

loaded = load('outputs/outputModelEPooled.mat');
data = struct();
data.t = loaded.tData;
data.r = loaded.rData;
data.gelPerc = [0,0.5,1];

loaded = load('outputs/outputModelA.mat');
fitted = struct();
fitted.t = loaded.reducedOutput.ts;
fitted.r = loaded.reducedOutput.rs;

save('outputs/ModelADataForFig.mat','data','fitted')