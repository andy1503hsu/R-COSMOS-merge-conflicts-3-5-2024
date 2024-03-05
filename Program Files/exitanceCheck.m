function model = exitanceCheck(model)
%% Exitance
    model.exitance.save = false; % save thermal/solar exitance results
    if model.exitance.save
    
    % Hemispherical angular grid for trajectory binning
    model.exitance.numZenithBins = 6; % Encompassing [0, pi/2)
    model.exitance.numAzimuthBins = 24; % Encompassing [0, 2pi)

    % Wavelength bounds for exitance detection
    model.exitance.solarBounds = [0 5]; % [microns]
    model.exitance.thermalBounds = [0 300]; % [microns]

    % Number of uniformally spaced bins within wavelength bounds
    model.exitance.numSolarBins = 50;
    model.exitance.numThermalBins = 225;

    % Exitance resolution (typically > model.resolution values)
    model.exitance.sourcePhotons = 2e6; % Light source photons
    model.exitance.thermalPhotons = 2e5; % Thermal photons per facet
    
    end
    
end