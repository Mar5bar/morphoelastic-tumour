function [vals, params] = fitFuncModelE(vecParams, tData, params)

    % Set values of params
    params.kappa = vecParams(1);

    % Run the simulation quietly.
    output = runSim(params, false);

    % Evaluate radii at the requested times.
    vals = interp1(output.ts, output.rs(:,end), tData);

    % Convert to um from m.
    vals = vals * 1e6;

end