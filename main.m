% Set up params in a struct.
params = struct(...  
        'T', 30,... % Final time.
        'nR', 1000,... % Spatial discretisation.
        'nT', 5000,... % Time discretisation
        'initialRadius', 1, ... % Initial dimensional radius of the spheroid.
        'radialStressBoundaryConstant', 0,... % Stress boundary condition constant `kappa'.
        'nutrientConsumptionRate', 0.15,... % Rate of nutrient consumption by tissue.
        'necrosisThreshold', 0.8,... % Necrosis threshold.
        'radialStressIntegrandThreshold', 0.05,... % Radial threshold for using Taylor expansion of radial stress integrand.
        'elasticStretchIntegrandThreshold', 0.05,... % Radial threshold for using Taylor expansion of elastic stretch.
        'stressGrowthThreshold', -1,... % Threshold below which growth is stopped by stress
        'stressResponse', 'global - boundary',... % Type of stress response. One of 'local', 'global - boundary', 'global - min', or 'none'.
        'stressType', 'radial',... % Type of stress to use in nAux. One of 'radial', 'hoop', or 'bulk'.
        'nAuxSelector', 2, ... % Selects the approxpriate auxiliary stress modifier. One of 1 or 2.
        'growthLawSelector', 2, ... % Selects the approxpriate growth law. One of 1 or 2.
        'necroticDecay', false ... % Does the necrotic material decay away?
);
% Compute the threshold radius of the spheroid, after which a core of zero
% nutrient forms.
params.outerBoundaryNutrientThreshold = sqrt(6/params.nutrientConsumptionRate);

% Initial discretisation.
RsMinusB = zeros(params.nT,params.nR); RsMinusB(1,:) = linspace(0,params.initialRadius,params.nR) - params.initialRadius;
ts = linspace(0,params.T,params.nT);

% Computed values to store.
rs = zeros(params.nT,params.nR);
growthStretches = zeros(params.nT,params.nR); growthStretches(1,:) = 1;
growthRates = zeros(params.nT,params.nR);
radii = zeros(params.nT,1);
radialStresses = zeros(params.nT,params.nR);
hoopStresses = zeros(params.nT,params.nR);
bulkStresses = zeros(params.nT,params.nR);
elasticStretches = zeros(params.nT,params.nR);
nutrients = zeros(params.nT,params.nR);
necroticRadii = zeros(params.nT,1)-Inf;

clear('textprogressbar.m')
textprogressbar('Simulating growth: ');
for tInd = 1 : params.nT
    textprogressbar(100*tInd/params.nT);

    % Compute the Eulerian radial coordinates.
    rs(tInd,:) = computeRadialCoord(RsMinusB(tInd,:),growthStretches(tInd,:),params);

    % Compute the radial, hoop, and bulk stresses in the current configuration.
    [radialStresses(tInd,:), hoopStresses(tInd,:), bulkStresses(tInd,:), elasticStretches(tInd,:)] = ...
            computeStresses(RsMinusB(tInd,:),rs(tInd,:),growthStretches(tInd,:),params);

    % Compute the current nutrient concentration.
    nutrients(tInd,:) = computeNutrient(rs(tInd,:),params);

    % Compute the current Eulerian radius of the necrotic core.
    necroticRadii(tInd) = computeNecroticRadius(necroticRadii(1:tInd-1),rs(tInd,:),nutrients(tInd,:),params);

    % Compute the growth rate.
    switch params.stressType
    case 'radial'
        growthRates(tInd,:) = computeGrowthRate(radialStresses(tInd,:),nutrients(tInd,:),rs(tInd,:),necroticRadii(tInd),params);
    case 'hoop'
        growthRates(tInd,:) = computeGrowthRate(hoopStresses(tInd,:),nutrients(tInd,:),rs(tInd,:),necroticRadii(tInd),params);
    case 'bulk'
        growthRates(tInd,:) = computeGrowthRate(bulkStresses(tInd,:),nutrients(tInd,:),rs(tInd,:),necroticRadii(tInd),params);
    end

    %% If there will be a next step
    if tInd < params.nT

        % Compute a new discretisation.
        RsMinusB(tInd+1,:) = interp1(rs(tInd,:), RsMinusB(tInd,:), linspace(rs(tInd,1),rs(tInd,end),params.nR));

        % Map current r, stress, nutrient, growthStretch and growthRate onto the new
        % discretisation, storing them temporarily.
        rRemapped = interp1(RsMinusB(tInd,:), rs(tInd,:), RsMinusB(tInd+1,:));
        stressRemapped = interp1(RsMinusB(tInd,:), radialStresses(tInd,:), RsMinusB(tInd+1,:));
        nutrientRemapped = interp1(RsMinusB(tInd,:), nutrients(tInd,:), RsMinusB(tInd+1,:));
        growthStretchRemapped = interp1(RsMinusB(tInd,:), growthStretches(tInd,:), RsMinusB(tInd+1,:));
        growthRateRemapped = interp1(RsMinusB(tInd,:), growthRates(tInd,:), RsMinusB(tInd+1,:));

        % Timestep gamma from tInd to tInd+1 using the mapped quantities from above.
        growthStretches(tInd+1,:) = growthStretchRemapped + ...
                                    (ts(tInd+1) - ts(tInd)) * growthStretchRemapped.*growthRateRemapped;
    end

end
textprogressbar('\nDone.');

% Postprocessing.
radii(:) = rs(:,end);

% Plotting.
plot_spheroid(ts, rs, growthStretches, growthRates, radialStresses, nutrients, params, 1);
plot_spheroid(ts, rs, growthStretches, growthRates, radialStresses, nutrients, params, params.nT);
plot_evolution(ts, rs, growthStretches, growthRates, radialStresses, nutrients, necroticRadii, params);


function r = computeRadialCoord(RMinusB,growthStretch,params)
%% Compute the Eulerian radial coordinates from the growth stretches at the
%% material points in RMinusB.
    r = (3*cumtrapz(RMinusB,growthStretch.^3.*(RMinusB+params.initialRadius).^2)).^(1/3);
end

function [radialStress, hoopStress, bulkStress, elasticStretch] = computeStresses(RMinusB,r,growthStretch,params)
%% Compute the radial and hoop stresses from the radial coordinates and growth
%% stretches at the material points in RMinusB.

    % Compute the radial stress first.
    integrand = 2*growthStretch.*(r.^6 - growthStretch.^6.*(RMinusB+params.initialRadius).^6)./r.^7;
    % Ignore the value at R=0, which is NaN.
    integrand(1) = 0;

    % The integrand is poorly computed around R = 0, so we use a two-term
    % Taylor expansion of the integrand around this point. We use the Taylor
    % expansion until we reach a threshold of R = params.radialStressIntegrandThreshold.
    growthStretchPrime = gradient(growthStretch, RMinusB);
    growthStretchPrimePrime = gradient(growthStretchPrime, RMinusB);
    approxIntegrand = 2*(-3/2 * growthStretchPrime(1)/gamma(1) + ...
                     (3/80 * (growthStretchPrime(1)/growthStretch(1))^2 - ...
                     6/5 * growthStretchPrimePrime(1)/growthStretch(1))*(RMinusB+params.initialRadius));
    integrandComposite = approxIntegrand.*(RMinusB<=(params.radialStressIntegrandThreshold-params.initialRadius)) + ...
                        integrand.*(RMinusB>(params.radialStressIntegrandThreshold-params.initialRadius));
    integral = cumtrapz(RMinusB,integrandComposite);

    radialStress = (integral - integral(end)) - params.radialStressBoundaryConstant*(r(end) - params.initialRadius)/params.initialRadius;

    % Compute the hoop stress from the radial stress.
    elasticStretch = r  ./ ((RMinusB+params.initialRadius).*growthStretch);
    % Ignore the value at R=0, which is NaN.
    elasticStretch(1) = 1;
    
    % As before, we'll use a two-term Taylor expansion to compute the
    % elastic stretch near R=0.
    approxElasticStretch = 1 - (growthStretchPrime(1)/gamma(1)).*(RMinusB+params.initialRadius)/4;
    elasticStretch = approxElasticStretch.*(RMinusB<=(params.elasticStretchIntegrandThreshold-params.initialRadius)) + ...
                        elasticStretch.*(RMinusB>(params.elasticStretchIntegrandThreshold-params.initialRadius));

    hoopStress = radialStress + (elasticStretch.^2 - 1./elasticStretch.^4);

    % Compute the bulk stress from the radial and hoop stress.
    bulkStress = radialStress + 2*hoopStress;

end

function nutrient = computeNutrient(r,params)
%% Compute the nutrient concentration at the material points of r.
    b = r(end);
    if b <= params.outerBoundaryNutrientThreshold
        % If no core of zero nutrient has developed:
        nutrient = 1 - params.nutrientConsumptionRate * (b^2 - r.^2) / 6;
    else
        % Else, a core of zero nutrient has developed: Find the root of the
        % polynomial that defines the radius of the core. We'll do this
        % crudely, then pass our guess into fzero.
        as = linspace(0,b,1e4);
        f = @(a) (2 * a.^3 - 3 * a.^2 * b + b^3) * params.nutrientConsumptionRate - 6 * b;
        fs = f(as);
        guessInd = find(fs(1:end-1).*fs(2:end)<=0,1,'first');
        % If no root has been found, the root must be between as(1) and as
        % (2). If so, we repeat the above calculation. The root is unlikely
        % to be within 1e-8 of zero, which is when this approach will fail.
        if isempty(guessInd)
            as = linspace(0,as(2),1e4);
            fs = f(as);
            guessInd = find(fs(1:end-1).*fs(2:end)<=0,1,'first');
        end
        rThreshold = fzero(f, as(guessInd:guessInd+1));
        nutrient = 0*r;
        mask = r > rThreshold;
        nutrient(mask) =params.nutrientConsumptionRate * r(mask).^2 / 6 + ...
                        (params.nutrientConsumptionRate * rThreshold^3 / 3)./r(mask) - ...
                        params.nutrientConsumptionRate * (2 * rThreshold^3 + b^3) / (6 * b) + ...
                        1;
    end
end


function necroticRadius = computeNecroticRadius(previousRadii,r,nutrient,params)
%% Returns the current radius of the necrotic core.
    % Compute the nutrient-determined candidate radius.
    nutrientDeterminedRadius = max([r(find(nutrient<params.necrosisThreshold,1,'last')),-Inf]);
    if params.necroticDecay
        % If the necrotic tissue can decay, we define necrosis via a Greenspan-like threshold.
        necroticRadius = nutrientDeterminedRadius;
    else
        % If the necrotic tissue cannot decay, we define necrosis as the max
        % of the nutrientDeterminedRadius and the previous radii.
        necroticRadius = max([nutrientDeterminedRadius; previousRadii]);
    end
end


function growthRate = computeGrowthRate(stress,nutrient,r,necroticRadius,params)
%% Returns the growthRate as a function of stress and nutrient.
    switch params.growthLawSelector
    case 1
        growthRate = growthStressResponse(stress,params) .* (nutrient - params.necrosisThreshold);
    case 2
        growthRateOne = growthStressResponse(stress,params) .* (nutrient - params.necrosisThreshold);
        if params.necroticDecay
            growthRateTwo = nutrient - params.necrosisThreshold;
        else
            growthRateTwo = 0;
        end
        mask = r > necroticRadius;
        growthRate = growthRateOne.*mask + growthRateTwo.*(1-mask);
    end
end

function n = growthStressResponse(stress,params)
%% Returns the growthStressResponse `n' as a function of stress.
    switch params.stressResponse
    case 'local'
        % For local stress response.
        n = nAux(stress,params);
    case 'global - boundary'
        % For global stress response to stress at the boundary:
        n = nAux(stress(end),params) * ones(size(stress));
    case 'global - min'
        % For global stress response to largest compressive stress:
        n = nAux(min(stress),params) * ones(size(stress));
    case 'none'
        % For no stress response.
        n = ones(size(stress));
    end
end

function nAux = nAux(stress,params)
%% Auxiliary function for the stress response.

    nAux = 0.*(stress < params.stressGrowthThreshold) + ...
         (1 - stress/params.stressGrowthThreshold).*(stress >= params.stressGrowthThreshold).*(stress < 0) + ...
         1.*(stress >= 0);
 end