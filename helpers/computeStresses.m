function [radialStress, hoopStress, bulkStress, elasticStretch] = computeStresses(RMinusB,r,growthStretch,params)
%% Compute the radial and hoop stresses from the radial coordinates and growth
%% stretches at the material points in RMinusB.

    % Compute the radial stress first.
    integrand = 2*growthStretch.*(r.^6 - growthStretch.^6.*(RMinusB+params.B).^6)./r.^7;
    % Ignore the value at R=0, which is NaN.
    integrand(1) = 0;

    % The integrand is poorly computed around R = 0, so we use a two-term
    % Taylor expansion of the integrand around this point. We use the Taylor
    % expansion until we reach a threshold of R = params.radialStressIntegrandThreshold.
    growthStretchPrime = gradient(growthStretch, RMinusB);
    growthStretchPrimePrime = gradient(growthStretchPrime, RMinusB);
    approxIntegrand = 2*(-3/2 * growthStretchPrime(1)/growthStretch(1) + ...
                     (3/80 * (growthStretchPrime(1)/growthStretch(1))^2 - ...
                     6/5 * growthStretchPrimePrime(1)/growthStretch(1))*(RMinusB+params.B));
    integrandComposite = approxIntegrand.*(RMinusB<=(params.radialStressIntegrandThreshold-params.B)) + ...
                        integrand.*(RMinusB>(params.radialStressIntegrandThreshold-params.B));
    integral = cumtrapz(RMinusB,integrandComposite);

    radialStress = params.mu*(integral - integral(end)) - params.kappa*(r(end) - params.B)/params.B;

    % Compute the hoop stress from the radial stress.
    elasticStretch = r  ./ ((RMinusB+params.B).*growthStretch);
    % Ignore the value at R=0, which is NaN.
    elasticStretch(1) = 1;
    
    % As before, we'll use a two-term Taylor expansion to compute the
    % elastic stretch near R=0.
    approxElasticStretch = 1 - (growthStretchPrime(1)/growthStretch(1)).*(RMinusB+params.B)/4;
    elasticStretch = approxElasticStretch.*(RMinusB<=(params.elasticStretchIntegrandThreshold-params.B)) + ...
                        elasticStretch.*(RMinusB>(params.elasticStretchIntegrandThreshold-params.B));

    hoopStress = radialStress + params.mu*(elasticStretch.^2 - 1./elasticStretch.^4);

    % Compute the bulk stress from the radial and hoop stress.
    bulkStress = radialStress + 2*hoopStress;

end