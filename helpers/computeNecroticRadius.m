function necroticRadius = computeNecroticRadius(previousRadii,r,nutrient,params)
%% Returns the current radius of the necrotic core.
    % Compute the nutrient-determined candidate radius.
    nutrientDeterminedRadius = max([r(find(nutrient<params.cHat,1,'last')),-Inf]);
    if contains('AB',params.model)
        % If the necrotic tissue can decay, we define necrosis via a Greenspan-like threshold.
        necroticRadius = nutrientDeterminedRadius;
    else
        % If the necrotic tissue cannot decay, we define necrosis as the max
        % of the nutrientDeterminedRadius and the previous radii.
        necroticRadius = max([nutrientDeterminedRadius; previousRadii]);
    end
end