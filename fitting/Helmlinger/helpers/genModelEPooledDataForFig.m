clear all

loaded = load('outputs/outputModelEPooled.mat');
data = struct();
data.t = loaded.tData;
data.r = loaded.rData;
data.gelPerc = [0,0.5,1];

fitted = struct();
fitted.t = loaded.tDataHighRes;
fitted.r = loaded.rsCell;

save('outputs/ModelEPooledDataForFig.mat','data','fitted')