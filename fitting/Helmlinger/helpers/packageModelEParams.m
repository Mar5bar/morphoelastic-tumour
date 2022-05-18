loaded = load('outputs/outputModelEPooled.mat');
fnames = {'0%','0.5%','1%'};
params = loaded.params;
params.k = loaded.x(1);
params.cHat = loaded.x(2);
params.lambda = loaded.x(3);
params.sigmaHat = loaded.x(4);
for i = 1 : 3

    % First dataset will always have kappa = 0.
    if i == 1
        params.kappa = 0;
    else
        params.kappa = loaded.x(3+i);
    end

    params.B = loaded.x(6+i);

    save(['fittedParameters/ModelEFittedParameters-',fnames{i},'.mat'],'params')

end