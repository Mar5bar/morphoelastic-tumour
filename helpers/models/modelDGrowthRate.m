function growthRate = modelDGrowthRate(nutrient,r,necroticRadius,radialStress,hoopStress,params)
%% Returns the growthRate 1/gamma * d(gamma)/dt corresponding to model D,
% recovering a robust steady state.

% In this model, growthRate = {k * (c - cHat) * n(min_R(sigma_r,sigma_theta)) if r >= necroticRadius,
%                             {0                                              if r <  necroticRadius.
globalMinStress = min([radialStress(:); hoopStress(:)]);
growthRateGrowing = params.k * (nutrient - params.cHat) .* nFun(globalMinStress,params);
growthRateNecrotic = 0;

% Identify necrotic tissue.
growingMask = r >= necroticRadius;

growthRate = growthRateGrowing.*growingMask + growthRateNecrotic.*(1-growingMask);

end