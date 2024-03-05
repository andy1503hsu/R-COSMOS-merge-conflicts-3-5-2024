function temp = binAbsorptionEvents(temp, aData, energyPerFacet)

    temp.thermalEmission = energyPerFacet; % HT boundary condition

    if temp.iter == 1  % First local day
        
        nodes = temp.meshNodes;
        elements = temp.meshElements;
        
        % Number of workers in parallel pool (used in binning algorithm)
        numWorkers = gcp().NumWorkers;

        if size(aData.solar, 1) == 0  % No solar absorption events
            temp.meshLightAbsorption = zeros(size(elements, 1), 1);
            temp.prevSolarAbsorption(:, temp.dayStep) = zeros(size(elements, 1), 1);
        else
            aEvents3D = padAbsorptionEvents(aData.solar, numWorkers);
            energyInMesh = binParallelized(aEvents3D, nodes, elements);

            temp.meshLightAbsorption = energyInMesh;
            temp.prevSolarAbsorption(:, temp.dayStep) = energyInMesh;
        end

        if temp.dayStep == 1 % Thermal radiation was also simulated
            energyInMesh = binParallelized(aData.thermal, nodes, elements);
            temp.thermalDistribution = energyInMesh ./ energyPerFacet';
            temp.meshThermalAbsorption = sum(energyInMesh, 2);
        else % Use percentage distribution to get thermal absorption
            fprintf("\n%30s %11.4f J\n", "Thermal energy emitted:", sum(energyPerFacet))
            energyInMesh = sum(temp.thermalDistribution .* energyPerFacet', 2); 
            fprintf("%30s %11.4f J %-20s\n", "Thermal energy absorbed:", sum(energyInMesh), "(from percentage distribution)")
            temp.meshThermalAbsorption = energyInMesh;
        end
        
    
    else % Not the first day, results will be reused

        temp.meshLightAbsorption = temp.prevSolarAbsorption(:, temp.dayStep);
        fprintf("\n%30s %11.4f J %-20s \n", "Solar energy absorbed:", sum(temp.meshLightAbsorption), "(from day 1)")
        
        fprintf("%30s %11.4f J\n", "Thermal energy emitted:", sum(energyPerFacet))
        energyInMesh = sum(temp.thermalDistribution .* energyPerFacet', 2); 
        fprintf("%30s %11.4f J %-20s\n", "Thermal energy absorbed:", sum(energyInMesh), "(from percentage distribution)")         
        temp.meshThermalAbsorption = energyInMesh;
    
    end
    %{
    testSolar = binningQuadTree(absorptionData.solar, nodes, elements, 1e7);
    testSnow = binningQuadTree(absorptionData.snow, nodes, elements, 1e7);
    %}
end