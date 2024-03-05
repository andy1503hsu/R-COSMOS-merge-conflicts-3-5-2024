function mieValuesCluster = mieCalculator(PMC)
    
    radius = PMC.grainRadii; % [microns]
    fprintf("Calculating single-scattering coefficients for r = %.2f microns:\n\n", radius);
    tic;

    % Complex refractive index of water-ice, from Warren and Brandt 2008
    wb2008 = readmatrix("WarrenBrandt2008_dataset.txt");
    % 1st column -- Wavelength, in microns
    % 2nd column -- Real component of refractive index
    % 3rd column -- Complex component of refractive index

    wb_wavelengths = wb2008(:, 1) * 1e-6; % [m]
    wb_refractIndex = wb2008(:, 2) + wb2008(:, 3)*1i;  % [complex]

    radius = radius * 1e-6; % [m] Convert radius from microns to meter
    
    % Get Mie values (Qext, omega, and g) for a *single* isolated water-ice sphere
    mieValuesSingle = getMieValues(wb_wavelengths, wb_refractIndex, radius);

    % Get Mie values for a *cluster* of water-ice spheres using the
    % delta-Eddington approximation, see Wiscombe and Warren 1980
    mieValuesCluster = deltaEddington(mieValuesSingle, PMC);

    fprintf("\nCompleted in %f s\n",toc);
end