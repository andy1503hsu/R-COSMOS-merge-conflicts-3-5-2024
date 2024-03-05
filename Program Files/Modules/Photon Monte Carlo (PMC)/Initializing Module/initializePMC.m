% Initialize/Preallocate values for PMC runs during a new global timestep.
% The snow surface remains constant during a global timestep, allowing the
% mesh nodes and elements to be cached. Furthermore, some matrices that
% are used to reduce Monte Carlo time are preallocated.

function temp = initializePMC(temp, sData, PMC)
    
    % Reset these values during first timestep of a new global step
    if temp.iter == 1 && temp.dayStep == 1
        %{
        % Get mesh nodes and elements
        if inPMC.type == "Planetary"
            msh = getPlanetaryMesh(sData, inPMC);
        elseif inPMC.type == "Ark Chamber"
            msh = getArkMesh(sData, inPMC);
            
        end
        temp.meshNodes = msh.Points;
        temp.meshElements = msh.ConnectivityList;
        %}
        
        % Get # of mesh elements and number of timesteps
        numElements = length(temp.meshElements);
        numDaySteps = size(temp.prevSolarAbsorption, 2);

        % Reuse solar absorption data from 1st day of global timestep
        temp.prevSolarAbsorption = zeros(numElements, numDaySteps);

        % Reuse thermal absorption data from 1st day of global timestep
        temp.thermalDistribution = zeros(numElements, length(temp.xS)-1);

        % Preallocate arrays for absorption
        temp.meshLightAbsorption = zeros(numElements, 1);
        temp.meshThermalAbsorption = zeros(numElements, 1);

        % Update solar irradiance values based on distance
        temp.specIrrad = specIrrad(PMC.lightSource, temp.lightDistance);
    end

end