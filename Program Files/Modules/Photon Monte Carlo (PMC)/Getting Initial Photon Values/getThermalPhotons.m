% Wall functionality is untested

function [thermalPhotons, energyPerFacet] = getThermalPhotons(temp, sData, PMC, type)
    
    % Physical constants
    c = 2.99792458e8; % [m/s] speed of light
    h = 6.62607015e-34; % [J-s] Plank's constant
    k = 1.380649e-23; % [J/K] Boltzmann's constant
    
    wavelengths = linspace(PMC.wavelengthBounds(1), ...
                           PMC.wavelengthBounds(2), 1e4) * 1e-6;

    if type == "Snow"
        
        bundlesPerSurface = PMC.thermalPhotons;
        totalBundles = bundlesPerSurface*length(sData.areaFacet);
    
    elseif type == "Wall" % SLI covering emits no photons
        
        bundlesPerSurface = PMC.wallPhotons;
        totalBundles = bundlesPerSurface*(length(sData.areaFacet) - 1);
    
    end
    
    nextEmptyBundle = 1;
    thermalPhotons = zeros(totalBundles, 16);
    energyPerFacet = zeros(length(sData.areaFacet), 1);
    
    % For every surface whose black body radiation we are concerned with...
    for index = 1:length(sData.areaFacet)
        if type == "Wall" && PMC.arkWalls.SLI_Index == index % Skip SLI covering
            continue
        end
        
        temperature = temp.tempProfile(index);        
        endIndex = nextEmptyBundle+bundlesPerSurface-1;
        
        % Energy irradiance (W/m^2/m) for one specific temperature
        irrad = PMC.emiFunc .* 2*pi*h*(c^2)./(wavelengths.^5)./...
                     (exp(h*c./wavelengths/k./temperature)-1); % Planck's law for black body radiation
        
        tInt = temp.time;
        photonsPerWavelength = irrad*tInt*sData.areaFacet(index).*wavelengths/h/c;
        
        % Sample and get wavelengths of each bundle, along with how many photons each bundle represents
        [bundleWavelengths, photonsPerBundle] = sampleDistribution(wavelengths, photonsPerWavelength, PMC.wavelengthBounds*1e-6, bundlesPerSurface);
        
        % Wavelengths and energies of simulated bundles
        thermalPhotons(nextEmptyBundle:endIndex, 8) = bundleWavelengths;
        thermalPhotons(nextEmptyBundle:endIndex, 9) = h*c./bundleWavelengths.*photonsPerBundle;
        energyPerFacet(index) = sum(h*c./bundleWavelengths.*photonsPerBundle);
        
        % Sample trajectory vectors from a Lambertian distribution
        zen = asin(sqrt(rand(bundlesPerSurface, 1))); % [rad] zenith angle
        azi = 2*pi*rand(bundlesPerSurface, 1); % [rad] azimuth angle
        
        if type == "Snow" % We want the counterclockwise-normal unit vector
            unit_XComp = -(sData.zf(index) - sData.zi(index)) / sData.areaFacet(index);
            unit_ZComp = (sData.xf(index) - sData.xi(index)) / sData.areaFacet(index);
        elseif type == "Wall" % Clockwise-normal unit vector
            unit_XComp = (sData.zf(index) - sData.zi(index)) / sData.areaFacet(index);
            unit_ZComp = -(sData.xf(index) - sData.xi(index)) / sData.areaFacet(index);
        end
        
        % Convert trajectory components from surface-relative components to "absolute" components
        if abs(unit_ZComp) > 0.9999
            thermalPhotons(nextEmptyBundle:endIndex, 5) = sin(zen).*cos(azi); % [-]
            thermalPhotons(nextEmptyBundle:endIndex, 6) = sin(zen).*sin(azi); % [-]
            thermalPhotons(nextEmptyBundle:endIndex, 7) = unit_ZComp*cos(zen)/abs(unit_ZComp); % [-]
        else
            thermalPhotons(nextEmptyBundle:endIndex, 5) = sin(zen)/sqrt(1-unit_ZComp^2).*(unit_XComp*unit_ZComp*cos(azi)) +...
                                                     unit_XComp*cos(zen); % [-]
            thermalPhotons(nextEmptyBundle:endIndex, 6) = sin(zen)/sqrt(1-unit_ZComp^2).*(unit_XComp*sin(azi)); % [-]
            thermalPhotons(nextEmptyBundle:endIndex, 7) = -sin(zen).*cos(azi)*sqrt(1-unit_ZComp^2) +...
                                                     unit_ZComp*cos(zen); % [-]
        end
    
        % Starting positions and regions
        randNums = rand(bundlesPerSurface, 1);
        
        % x-coordinates randomly sampled along surface
        thermalPhotons(nextEmptyBundle:endIndex, 1) = sData.xi(index) + (sData.xf(index) - sData.xi(index))*randNums;
        % y-coordinates uniformly and randomly sampled along unit depth
        thermalPhotons(nextEmptyBundle:endIndex, 2) = randNums;
        % z-coordinates randomly sampled along surface
        thermalPhotons(nextEmptyBundle:endIndex, 3) = sData.zi(index) + (sData.zf(index) - sData.zi(index))*randNums;
        
        if type == "Wall" && index == 1 % Photons that are emitted from the Ark Chamber floor directly under the snow start in the snow region
            thermalPhotons(nextEmptyBundle:endIndex, 4) = 2;
        else % Otherwise, photons start in ray-tracing region
            thermalPhotons(nextEmptyBundle:endIndex, 4) = 1;
        end
        
        nextEmptyBundle = endIndex + 1;
    end

    if type == "Snow" % Snow photons start on a surface
        thermalPhotons(:, 15) = 1;
    end

    %  Single-scattering coefficients and energy-related values
    thermalPhotons = addRemainingProperties(thermalPhotons, PMC);
end