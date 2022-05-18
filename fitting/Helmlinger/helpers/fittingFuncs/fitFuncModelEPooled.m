function [res, resCell] = fitFuncModelEPooled(vecParams, tData, rData, params, parWorkers)

    % Set values of params
    params.k = vecParams(1);
    params.cHat = vecParams(2);
    params.lambda = vecParams(3);
    params.sigmaHat = vecParams(4);

    res = zeros(sum(cellfun(@numel, tData)),1);
    resCell = cell(length(rData),1);
    ind = 1;

    if nargin < 5
        parWorkers = Inf;
    end

    parfor (i = 1 : length(rData), parWorkers)

        curParams = params;

        % First dataset will always have kappa = 0.
        if i == 1
            curParams.kappa = 0;
        else
            curParams.kappa = vecParams(3+i);
        end

        curParams.B = vecParams(6+i);

        % Run the simulation quietly.
        output = runSim(curParams, false);

        % Evaluate radii at the requested times.
        vals = interp1(output.ts, output.rs(:,end), tData{i});

        % Convert to um from m.
        vals = vals * 1e6;

        % Store the residual of vals vs rData{i}.
        resCell{i} = vals(:) - rData{i}(:);

    end

    res = cell2mat(resCell);

end