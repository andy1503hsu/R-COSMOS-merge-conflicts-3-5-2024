function spectralRadiance = binExitData(temp, eData, PMC, type)

    nZen = PMC.exitance.numZenithBins;
    nAzi = PMC.exitance.numAzimuthBins;

    if type == "thermal"
        nThermal = PMC.exitance.numThermalBins;
        bounds = PMC.exitance.thermalBounds*1e-6; % [m]
    
        spectralRadiance = zeros(nThermal, nZen, nAzi);
        wavelengthEdges = linspace(bounds(1), bounds(2), nThermal+1);
    
    elseif type == "solar"
        nSolar = PMC.exitance.numSolarBins;
        bounds = PMC.exitance.solarBounds*1e-6; % [m]

        spectralRadiance = zeros(nSolar, nZen, nAzi);
        wavelengthEdges = linspace(bounds(1), bounds(2), nSolar+1);
    end

    % Eliminate any exit events that aren't within wavelength bounds
    eData(eData(:, 1) >= wavelengthEdges(end) | eData(:, 1) < wavelengthEdges(1), :) = [];

    % Get zenith/azimuth edges and indices for all exit events
    [zenEdges, aziEdges, zenIds, aziIds] = binZenAzi(eData, nZen, nAzi);

    for i = 1:nZen
        for j = 1:nAzi
            relevantExits = eData(zenIds == i & aziIds == j, :);
            if isempty(relevantExits) % No photons within this angular bin
                spectralRadiance(:, i, j) = 0;
            else
                sr_ = getSpectralRadiance(relevantExits, wavelengthEdges, zenEdges, aziEdges, i, j, temp);
                spectralRadiance(:, i, j) = sr_;
            end
        end
    end

end