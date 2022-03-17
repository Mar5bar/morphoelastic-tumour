function output = runSim(params)

    % Initial discretisation.
    RsMinusB = zeros(params.nT,params.nR); RsMinusB(1,:) = linspace(0,params.B,params.nR) - params.B;
    ts = linspace(0,params.tFinal,params.nT);

    % Computed values to store.
    rs = zeros(params.nT,params.nR);
    growthStretches = zeros(params.nT,params.nR); growthStretches(1,:) = 1;
    growthRates = zeros(params.nT,params.nR);
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
        switch params.model
        case 'A'
            growthRates(tInd,:) = modelAGrowthRate(nutrients(tInd,:),params);
        case 'B'
            growthRates(tInd,:) = modelBGrowthRate(nutrients(tInd,:),rs(tInd,:),necroticRadii(tInd),radialStresses(tInd,:),params);
        case 'C'
            growthRates(tInd,:) = modelCGrowthRate(nutrients(tInd,:),rs(tInd,:),necroticRadii(tInd),radialStresses(tInd,:),params);
        case 'D'
            growthRates(tInd,:) = modelDGrowthRate(nutrients(tInd,:),rs(tInd,:),necroticRadii(tInd),radialStresses(tInd,:),hoopStresses(tInd,:),params);
        case 'E'
            growthRates(tInd,:) = modelEGrowthRate(nutrients(tInd,:),rs(tInd,:),necroticRadii(tInd),radialStresses(tInd,:),hoopStresses(tInd,:),params);
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

    % Package output into a structure.
    output = struct();
    output.ts = ts;
    output.RsMinusB = RsMinusB;
    output.rs = rs;
    output.growthStretches = growthStretches;
    output.growthRates = growthRates;
    output.radialStresses = radialStresses;
    output.hoopStresses = hoopStresses;
    output.bulkStresses = bulkStresses;
    output.elasticStretches = elasticStretches;
    output.nutrients = nutrients;
    output.necroticRadii = necroticRadii;
    output.params = params;

end