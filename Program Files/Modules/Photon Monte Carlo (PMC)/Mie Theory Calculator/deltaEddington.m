function newValues = deltaEddington(mieValues, PMC)
    % Converts the Mie values (extinction effiency Qext, single-scattering
    % albedo omega, and asymmetry parameter g) of a isolated water-ice
    % sphere to the values for a cluster of water-ice spheres. These are:
    %   Attenuation coefficient mu  [1/m]
    %   Single-scattering albedo omega [-]
    %   Asymmetry Parameter g [-]

    % This conversion is known as the delta-Eddington approximation. See
    % Wiscombe and Warren (1980) for more information.

    fprintf("%5s Running delta-Eddington approximation... ", "")
    wavelengths = mieValues(:, 1);
    Qext = mieValues(:, 2);
    omega = mieValues(:, 3);
    g = mieValues(:, 4);

    rhoIce = 917; % kg/m^3, density of water-ice at 273 K
    radius = PMC.grainRadii * 1e-6; % [m] Convert from microns to m

    % delta-Eddington approximation
    mu = 0.75 * PMC.density/rhoIce * Qext / radius .* (1 - omega.*g.^2);
    omega2 = ((1 - g.^2).*omega) ./ (1 - g.^2.*omega);
    g2 = g ./ (1 + g);

    newValues = [wavelengths mu omega2 g2];

    fprintf("done\n")
end