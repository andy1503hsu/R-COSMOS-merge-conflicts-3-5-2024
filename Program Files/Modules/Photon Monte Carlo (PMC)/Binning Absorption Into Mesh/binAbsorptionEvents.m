function binnedData = binAbsorptionEvents(temp, aData, energyPerFacet, cP)

    if temp.iter == 1  % First local day
        
        nodes = temp.meshNodes;
        elements = temp.meshElements;
        
        % Number of workers in parallel pool (used in binning algorithm)
        numWorkers = gcp().NumWorkers;

        if size(aData.solar, 1) == 0  % No solar absorption events
            binnedData.meshLightAbsorption = zeros(size(elements, 1), 1);
        else
            aEvents3D = padAbsorptionEvents(aData.solar, numWorkers);
            energyInMesh = binParallelized(aEvents3D, nodes, elements);
            binnedData.meshLightAbsorption = energyInMesh;
        end

        if temp.dayStep == 1 % Thermal radiation was also simulated
            energyInMesh = binParallelized(aData.thermal, nodes, elements);
            binnedData.thermalDistribution = energyInMesh ./ energyPerFacet';
            binnedData.meshThermalAbsorption = sum(energyInMesh, 2);
        elseif cP == 1 % Use percentage distribution to get thermal absorption
            fprintf("\n%30s %11.4f J\n", "Thermal energy emitted:", sum(energyPerFacet))
            energyInMesh = sum(temp.thermalDistribution .* energyPerFacet', 2); 
            fprintf("%30s %11.4f J %-20s\n", "Thermal energy absorbed:", sum(energyInMesh), "(from percentage distribution)")
            binnedData.meshThermalAbsorption = energyInMesh;
        end
        
    
    else % Not the first day, results will be reused

        binnedData.meshLightAbsorption = temp.prevSolarAbsorption(:, temp.dayStep);
        fprintf("\n%30s %11.4f J %-20s \n", "Solar energy absorbed:", sum(temp.meshLightAbsorption), "(from day 1)")
        
        fprintf("%30s %11.4f J\n", "Thermal energy emitted:", sum(energyPerFacet))
        energyInMesh = sum(temp.thermalDistribution .* energyPerFacet', 2); 
        fprintf("%30s %11.4f J %-20s\n", "Thermal energy absorbed:", sum(energyInMesh), "(from percentage distribution)")         
        binnedData.meshThermalAbsorption = energyInMesh;
    
    end

end