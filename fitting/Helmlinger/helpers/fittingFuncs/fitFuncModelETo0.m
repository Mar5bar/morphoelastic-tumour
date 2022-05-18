function [vals, params] = fitFuncModelE(vecParams, tData, params)

    % Set values of params
    params.k = vecParams(1);
    params.cHat = vecParams(2);
    params.lambda = vecParams(3);
    params.B = vecParams(4);
    params.sigmaHat = vecParams(5);

    % Recompute the dependent parts of params.
    params.L = sqrt(params.D * params.cInf/params.lambda); % The diffusive lengthscale.
    params.bHat = sqrt(6)*params.L; % The threshold radius of the tumour for full/partial perfusion.

    % Run the simulation quietly.
    output = runSim(params, false);

    % Evaluate radii at the requested times.
    vals = interp1(output.ts, output.rs(:,end), tData);

    % Convert to um from m.
    vals = vals * 1e6;

end