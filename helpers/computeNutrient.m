function nutrient = computeNutrient(r,params)
%% Compute the nutrient concentration at the material points of r.
    b = r(end);
    if b <= params.bHat
        % If no core of zero nutrient has developed:
        nutrient = params.cInf - params.lambda * (b^2 - r.^2) / (6 * params.D);
    else
        % Else, a core of zero nutrient has developed: Find the root of the
        % polynomial that defines the radius of the core. We'll do this
        % crudely, then pass our guess into fzero.
        as = linspace(0,b,1e4);
        f = @(a) (2 * a.^3 - 3 * a.^2 * b + b^3) * params.lambda - 6 * b * params.D * params.cInf;
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
        % At this point, we just set rThreshold to zero.
        if isempty(guessInd)
            rThreshold = 0;
        else
            rThreshold = fzero(f, as(guessInd:guessInd+1));
        end
        nutrient = 0*r;
        mask = r > rThreshold;
        nutrient(mask) = params.lambda * r(mask).^2 / (6 * params.D) + ...
                        (params.lambda * rThreshold^3 / (3 * params.D)) ./ r(mask) - ...
                        params.lambda * (2 * rThreshold^3 + b^3) / (6 * b * params.D) + ...
                        params.cInf;
    end
end