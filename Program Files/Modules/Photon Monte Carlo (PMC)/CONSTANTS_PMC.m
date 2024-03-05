% Many parameters for the PMC modules are more or less constants -- for
% most simulations, these values are never adjusted. As a result, to
% declutter mainUI.m, it was decided to put all of these psuedo-constants
% in a separate file.

% Of course, these values can still be changed if need be, but the purpose 
% of decluttering mainUI.m still holds.

function PMC = CONSTANTS_PMC(PMC) %#codegen

    % Planetary case only parameters
    if PMC.type == "Planetary"
        % Are there any?
    % Ark Chamber case only parameters
    elseif PMC.type == "Ark Chamber"
         % [-] wall reflectivity + wall emissivity must = 1
        PMC.arkWalls.reflectivityFunction = @(wavelength) wallReflectivityArk(wavelength);
        PMC.arkWalls.emissivityFunction = @(wavelength) 1 - PMC.arkWalls.reflectivityFunction(wavelength); % [-]
        PMC.arkWalls.SLIReflectivity = 0.95; % Berisford 2018, page 5
        
        PMC.wallPhotons = 1e5; % [photons] Photons per wall segment
        PMC.arkDimensions = [1.2, 0.44]; % [m] Width, height
        PMC.SLIWidth = 0.279; % [m] SLI width on Ark ceiling
        PMC.LED_Dimensions = [0.011 0.136]; %[m] Width, height from ceiling
        
        PMC.arkWalls.wallsPerSide = 10;
        PMC.arkWalls.wallsFromCornerToSLI = 10;
        PMC.arkWalls.floorTemperature = 157; % [K]
        PMC.arkWalls.SLITemperature = 225; % [K]
        PMC.arkWalls.tempGradientType = "linear";

        [wallTemperatures, wallCoords, SLI_index] = arkWallsGradient(PMC.arkDimensions, PMC.arkWalls, PMC.SLIWidth);
        
        PMC.arkWalls.temperatures = wallTemperatures;
        PMC.arkWalls.xS = wallCoords.xS;
        PMC.arkWalls.zS = wallCoords.zS;
        PMC.arkWalls.SLI_Index = SLI_index;
    end
end