function growthRate = modelAGrowthRate(nutrient,params)
%% Returns the growthRate 1/gamma * d(gamma)/dt corresponding to model A, a
% minimal nutrient-limited growth model.

% In this model, growthRate = k * (c - cHat) for all tissue.
growthRate = params.k * (nutrient - params.cHat);

end