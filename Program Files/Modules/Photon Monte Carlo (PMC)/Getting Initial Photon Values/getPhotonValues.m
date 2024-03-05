function [photons, energyPerFacet] = getPhotonValues(temp, sData, PMC) %#codegen

    PMC.emiFunc = PMC.emissivity; % emissivity of snow

    % During first day of a global timestep...
    if temp.iter == 1

        % Get "non-changing" photons
            % For Planetary cases, these are solar photons
            % For Ark cases, these are LED and Ark Walls photons
        if PMC.type == "Planetary"
            solarPhotons = getSolarPhotons(temp, PMC);
        elseif PMC.type == "Ark Chamber"
            disp("ARK CHAMBER IN NEW PMC IS UNTESTED")
            solarPhotons = getLEDPhotons(temp, PMC);
            wallPhotons = getThermalPhotons(temp, temp.wall_sData, PMC, "Wall");
            solarPhotons = [solarPhotons; wallPhotons];
        end
    
        % Get thermal photons, only for first step
        if temp.dayStep == 1
            [thermalPhotons, ...
             energyPerFacet] = getThermalPhotons(temp, sData, PMC, "Snow");      
        else
            thermalPhotons = zeros(0, 16);
            energyPerFacet = getFacetEmission(temp, sData, PMC, "Snow");
        end

        photons.values = [solarPhotons; thermalPhotons];
        photons.thermalStart = size(solarPhotons, 1) + 1;
    
    else  % >= 2nd day of a global timestep
    
        % No photons (solar, thermal, wall) are traced
        energyPerFacet = getFacetEmission(temp, sData, PMC, "Snow");
        photons.values = zeros(0, 16);
        photons.thermalStart = NaN;
    
    end

end