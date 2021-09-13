% Set up params in a struct.
params = struct(...  
        'T', 30,... % Final time.
        'nR', 1000,... % Spatial discretisation.
        'nT', 5000,... % Time discretisation
        'stressBoundaryConstant', 0,... % Stress boundary condition constant `kappa'.
        'nutrientConsumptionRate', 0.15,... % Rate of nutrient consumption by tissue.
        'necrosisThreshold', 0.8,... % Necrosis threshold.
        'stressIntegrandThreshold', 1e-1,... % Radial threshold for using Taylor expansion of integrand.
        'stressGrowthThreshold', -1,... % Threshold below which growth is stopped by stress
        'globalStressResponse', true... % Boolean to signify global (true) or local (false) stress response.
);
% Compute the threshold radius of the spheroid, after which a core of zero
% nutrient forms.
params.outerBoundaryNutrientThreshold = sqrt(6/params.nutrientConsumptionRate);

% Initial discretisation.
RsMinusOne = zeros(params.nT,params.nR); RsMinusOne(1,:) = linspace(0,1,params.nR) - 1;
ts = linspace(0,params.T,params.nT);

% Computed values to store.
rs = zeros(params.nT,params.nR);
growthStretches = zeros(params.nT,params.nR); growthStretches(1,:) = 1;
growthRates = zeros(params.nT,params.nR);
radii = zeros(params.nT,1);
stresses = zeros(params.nT,params.nR);
nutrients = zeros(params.nT,params.nR);

textprogressbar('Simulating growth: ');
for tInd = 1 : params.nT
    textprogressbar(100*tInd/params.nT);

    % Compute the Eulerian radial coordinates.
    rs(tInd,:) = computeRadialCoord(RsMinusOne(tInd,:),growthStretches(tInd,:));

    % Compute the radial stresses in the current configuration.
    stresses(tInd,:) = computeStress(RsMinusOne(tInd,:),rs(tInd,:),growthStretches(tInd,:),params);

    % Compute the current nutrient concentration.
    nutrients(tInd,:) = computeNutrient(rs(tInd,:),params);

    % Compute the growth rate.
    growthRates(tInd,:) = computeGrowthRate(stresses(tInd,:),nutrients(tInd,:),params);

    %% If there will be a next step
    if tInd < params.nT

        % Compute a new discretisation.
        RsMinusOne(tInd+1,:) = interp1(rs(tInd,:), RsMinusOne(tInd,:), linspace(rs(tInd,1),rs(tInd,end),params.nR));

        % Map current r, stress, nutrient, growthStretch and growthRate onto the new
        % discretisation, storing them temporarily.
        rRemapped = interp1(RsMinusOne(tInd,:), rs(tInd,:), RsMinusOne(tInd+1,:));
        stressRemapped = interp1(RsMinusOne(tInd,:), stresses(tInd,:), RsMinusOne(tInd+1,:));
        nutrientRemapped = interp1(RsMinusOne(tInd,:), nutrients(tInd,:), RsMinusOne(tInd+1,:));
        growthStretchRemapped = interp1(RsMinusOne(tInd,:), growthStretches(tInd,:), RsMinusOne(tInd+1,:));
        growthRateRemapped = interp1(RsMinusOne(tInd,:), growthRates(tInd,:), RsMinusOne(tInd+1,:));

        % Timestep gamma from tInd to tInd+1 using the mapped quantities from above.
        growthStretches(tInd+1,:) = growthStretchRemapped + ...
                                    (ts(tInd+1) - ts(tInd)) * growthStretchRemapped.*growthRateRemapped;
    end

end
textprogressbar('\nDone.');

% Postprocessing.
radii(:) = rs(:,end);

% Plotting.
plot_spheroid(ts, rs, growthStretches, growthRates, stresses, nutrients, params, 1);
plot_spheroid(ts, rs, growthStretches, growthRates, stresses, nutrients, params, params.nT);
plot_evolution(ts, rs, growthStretches, growthRates, stresses, nutrients, params);


function r = computeRadialCoord(RMinusOne,growthStretch)
%% Compute the Eulerian radial coordinates from the growth stretches at the
%% material points in RMinusOne.
    r = (3*cumtrapz(RMinusOne,growthStretch.^3.*(RMinusOne+1).^2)).^(1/3);
end

function stress = computeStress(RMinusOne,r,growthStretch,params)
%% Compute the radial stresses from the radial coordinates and growth
%% stretches at the material points in RMinusOne.
    integrand = 2*growthStretch.*(r.^6 - growthStretch.^6.*(RMinusOne+1).^6)./r.^7;
    % Ignore the value at R=0, which is NaN.
    integrand(1) = 0;

    % The integrand is poorly computed around R = 0, so we use a two-term
    % Taylor expansion of the integrand around this point. We use the Taylor
    % expansion until we reach a threshold of R = params.stressIntegrandThreshold.
    growthStretchPrime = gradient(growthStretch, RMinusOne);
    growthStretchPrimePrime = gradient(growthStretchPrime, RMinusOne);
    approxIntegrand = 2*(-3/2 * growthStretchPrime(1)/gamma(1) + ...
                     (3/80 * (growthStretchPrime(1)/growthStretch(1))^2 - ...
                     6/5 * growthStretchPrimePrime(1)/growthStretch(1))*(RMinusOne+1));
    integrandComposite = approxIntegrand.*(RMinusOne<=(params.stressIntegrandThreshold-1)) + ...
                        integrand.*(RMinusOne>(params.stressIntegrandThreshold-1));
    integral = fliplr(cumtrapz(fliplr(RMinusOne),fliplr(integrandComposite)));

    stress = - params.stressBoundaryConstant*(r(end) - 1) - integral;
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


function growthRate = computeGrowthRate(stress,nutrient,params)
%% Returns the growthRate as a function of stress and nutrient.
    growthRate = growthStressResponse(stress,params) .* (nutrient - params.necrosisThreshold);
end

function n = growthStressResponse(stress,params)
%% Returns the growthStressResponse `n' as a function of stress.
    if params.globalStressResponse
        % For global stress response to stress at the boundary:
        n = nAux(stress(end),params) * ones(size(stress));
    else
        % For local stress response.
        n = nAux(stress,params);
    end
end

function nAux = nAux(stress,params)
%% Auxiliary function for the stress response.

    nAux = 0.*(stress < params.stressGrowthThreshold) + ...
         (1 - stress/params.stressGrowthThreshold).*(stress >= params.stressGrowthThreshold).*(stress < 0) + ...
         1.*(stress >= 0);
 end