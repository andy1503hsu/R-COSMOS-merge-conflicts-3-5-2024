% Wall functionality is untested

% This function is exactly identical to getThermalPhotons.m, but only the
% emitted energy per facet is outputted (no to-be-simulated photon
% bundles).

function energyPerFacet = getFacetEmission(temp, sData, PMC, type)
    
    % Physical constants
    c = 2.99792458e8; % [m/s] speed of light
    h = 6.62607015e-34; % [J-s] Plank's constant
    k = 1.380649e-23; % [J/K] Boltzmann's constant
    
    wavelengths = linspace(PMC.wavelengthBounds(1), ...
                           PMC.wavelengthBounds(2), 1e4) * 1e-6;

    if type == "Snow"
        
        bundlesPerSurface = PMC.thermalPhotons;
    
    elseif type == "Wall" % SLI covering emits no photons
        
        bundlesPerSurface = PMC.wallPhotons;
    
    end
    
    energyPerFacet = zeros(length(sData.areaFacet), 1);
    
    % For every surface whose radiation we are concerned with...
    for index = 1:length(sData.areaFacet)

        if type == "Wall" && PMC.arkWalls.SLI_Index == index % Skip SLI covering
            continue
        end
        
        temperature = temp.tempProfile(index);        
        
        % Energy irradiance (W/m^2/m) for one specific temperature
        irrad = PMC.emiFunc .* 2*pi*h*(c^2)./(wavelengths.^5)./...
                     (exp(h*c./wavelengths/k./temperature)-1); % Planck's law for black body radiation
        
        tInt = temp.time;
        photonsPerWavelength = irrad*tInt*sData.areaFacet(index).*wavelengths/h/c;
        
        % Sample and get wavelengths of each bundle, along with how many photons each bundle represents
        [bundleWavelengths, photonsPerBundle] = sampleDistribution(wavelengths, photonsPerWavelength, PMC.wavelengthBounds * 1e-6, bundlesPerSurface);
        
        % Total energy emitted from this facet
        energyPerFacet(index) = sum(h*c./bundleWavelengths.*photonsPerBundle);
        
    end
end