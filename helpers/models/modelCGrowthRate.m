function growthRate = modelCGrowthRate(nutrient,r,necroticRadius,radialStress,params)
%% Returns the growthRate 1/gamma * d(gamma)/dt corresponding to model C,
% a different approach to necrosis.

% In this model, growthRate = {k * (c - cHat) * n(sigma_r)   if r >= necroticRadius,
%                             {0                             if r <  necroticRadius.

growthRateGrowing = params.k * (nutrient - params.cHat) .* nFun(radialStress,params);
growthRateNecrotic = 0;

% Identify necrotic tissue.
growingMask = r >= necroticRadius;

growthRate = growthRateGrowing.*growingMask + growthRateNecrotic.*(1-growingMask);

end