function growthRate = modelBGrowthRate(nutrient,r,necroticRadius,radialStress,params)
%% Returns the growthRate 1/gamma * d(gamma)/dt corresponding to model B,
% coupling growth to stress.

% In this model, growthRate = {k * (c - cHat) * n(sigma_r)   if c >= cHat,
%                             {k * (c - cHat)                if c <  cHat.

growthRateGrowing = params.k * (nutrient - params.cHat) .* nFun(radialStress,params);
growthRateNecrotic = params.k * (nutrient - params.cHat);

% Identify necrotic tissue.
growingMask = r >= necroticRadius;

growthRate = growthRateGrowing.*growingMask + growthRateNecrotic.*(1-growingMask);

end