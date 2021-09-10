% Set up params in a struct.
params = struct(...  
        'T', 30,... % Final time.
        'nR', 1000,... % Spatial discretisation.
        'nT', 5000,... % Time discretisation
        'stresBoundaryConstant', 0,... % Stress boundary condition constant `kappa'.
        'nutrientConsumptionRate', 0.15,... % Rate of nutrient consumption by tissue.
        'necrosisThreshold', 0.8,... % Necrosis threshold.
        'stressIntegrandThreshold', 1e-1,... % Radial threshold for using Taylor expansion of integrand.
        'stressGrowthThreshold', -1,... % Threshold below which growth is stopped by stress
        'globalStressResponse', true... % Boolean to signify global (true) or local (false) stress response.
);

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
    nutrients(tInd,:) = 1 - params.nutrientConsumptionRate * (rs(tInd,end)^2 - rs(tInd,:).^2) / 6;

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

    stress = - params.stresBoundaryConstant*(r(end) - 1) - integral;
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