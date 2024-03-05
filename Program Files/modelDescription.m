function modelDescription(model)
%       Displays the model inputs: time discretization, light source,
% wavelength bounds, material properties, and numerical resolution.

    % Simulation type
    fprintf('Simulation type: %s\n',model.type)

    % Clear a log from storage and in-memory and start a new recording
    diary off
    if isfile('R-COSMOS.log')
        delete R-COSMOS.log % delete physical file
    end
    diary R-COSMOS.log % start a new log

    % Remove temporary files from a previous simuation
    if isfolder('Temporary')
        [~, ~] = rmdir('Temporary', 's'); % remove all folders and subfolders
    end

    % Create new temporary folders to store compilation files
    mkdir Temporary restartHead % restart head
    mkdir Temporary PMC % Photon Monte Carlo
    mkdir Temporary BMC % Blackbody Monte Carlo
    mkdir Temporary SMC % Sublimation Monte Carlo
    addpath(genpath("Temporary")) % add Temporary folder to the search path

    % Display user inputs
    % Time
    fprintf('\nTime\n')
    fprintf('  Geologic periods: %d\n', model.resolution.geologicSteps)
    fprintf('  Total simulation time: '); timeDisplay(model.time.total)
    fprintf('  Maximum iteration days: %d\n', model.resolution.maxDays)
    fprintf('  Day steps: %d\n', model.resolution.daySteps)
    fprintf('  Day length: '); timeDisplay(timePeriod(rotationPeriod(model.body.rotation)))
    fprintf('  Steady-state tolerance: %.2f %% \n', model.resolution.tolerance*100)

    % Orbit
    fprintf('\nOrbit\n')
    fprintf('  Sidereal period: '); timeDisplay(siderealPeriod(model.orbit.period))
    fprintf('  Orbit steps: %d\n', model.resolution.orbitSteps)
    fprintf('  Axial tilt: %.2f deg\n', axialTilt(model.orbit.axialTilt))
    fprintf('  Eccentricity: %.2f\n', eccentricity(model.orbit.eccentricity))
    fprintf('  Initial true anomaly: %d deg\n', model.orbit.initialTrueAnomaly)

    % Light parameters
    fprintf('\nLight\n')
    fprintf('  Light source: %s\n', model.light.source)
    fprintf('  Wavelength bounds: %.2f', model.light.wavelengthBounds(1))
    fprintf(' to %.0f microns\n', model.light.wavelengthBounds(2))

    % Surface
    fprintf('\nSurface\n')
    fprintf('  Type: %s\n', model.surface.type)
    fprintf('  Height: %.4f m\n', model.surface.height)
    fprintf('  Width: %.4f m\n', model.surface.width)
    fprintf('  Depth: %.4f m\n', model.surface.depth)
    fprintf('  Number: %d\n', model.surface.number)
    fprintf('  Facets: %d\n', model.surface.facets)
    fprintf('  Roughness: %.2f %%\n', model.surface.roughness)
    fprintf('  Modifier: %d\n', model.surface.modifier)
    fprintf('  Latitude: %.2f deg\n', model.surface.latitude)
    fprintf('  Orientation: %.2f deg West\n', surfaceOrientation(model.surface.orientation))

    % Material properties
    fprintf('\nMaterial properties\n')
    fprintf('  Species: %s\n', model.material.species)
    fprintf('  Density: %.2f kg/m^3\n', model.material.density)
    fprintf('  Grain radii: %d microns\n', model.material.grainRadii)
    fprintf('  Thermal inertia: %.2f J/m^2-K-s^0.5\n', model.material.thermalInertia)
    fprintf('  Initial temperature: %.2f K\n', model.material.initialTemperature)

    % Numerical resolution
    fprintf('\nNumerical resolution\n')
    fprintf('  Gradient growth: %.2f\n', model.resolution.Hgrad)
    fprintf('  Light photons: %.1g\n', model.resolution.sourcePhotons)
    fprintf('  Thermal photons per facet: %.1g\n', model.resolution.thermalPhotons)
    fprintf('  Molecules per facet: %.1g\n\n', model.resolution.molecules)

end
