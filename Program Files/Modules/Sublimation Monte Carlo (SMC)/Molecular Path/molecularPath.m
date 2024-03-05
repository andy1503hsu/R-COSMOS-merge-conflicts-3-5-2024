function location = molecularPath(sData,launchData,SMC) %#codegen
%       Simulates parallel ray-tracing according to the type of simulation
%   selected. The parallel portion of the simulation is included with each
%   module.
    
    % Simulation selector (boundary conditions included)
    % Real-world (meter-scale)
    if SMC.type == "Open"
        location = openSMC(sData,launchData);
        
    % Microscopic laboratory samples, ice crystals, comets
    elseif SMC.type == "Free"
        location = freeSMC(sData,launchData);
        
    % Planetary large-scale
    elseif SMC.type == "Planetary"
        if SMC.gravity % to conserve mass in the system
            location = planetaryWithGravitySMC(sData,launchData,SMC);
            
        else % for regions of net sublimation
            location = planetarySMC(sData,launchData);
        end
        
    % ARK chamber | Ocean Worlds Lab
    elseif SMC.type == "ArkChamber"
        location = ArkChamberSMC(sData,launchData);
        
    else
        error('Simulation type ''%s'' is not installed.',SMC.type)
    end
end