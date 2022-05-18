function r = computeRadialCoord(RMinusB,growthStretch,params)
%% Compute the Eulerian radial coordinates from the growth stretches at the
%% material points in RMinusB.
    r = (3*cumtrapz(RMinusB,growthStretch.^3.*(RMinusB+params.B).^2)).^(1/3);
end