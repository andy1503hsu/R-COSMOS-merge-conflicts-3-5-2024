function temp = radiationSolver(temp, sData, PMC)

    % Get values of to-be-simulated photons (if any), as well as
    % the energy emitted by each facet
    [photons, energyPerFacet] = getPhotonValues(temp, sData, PMC);
    temp.thermalEmission = energyPerFacet; % HT boundary condition

    if temp.iter == 1  % Partitioning occurs on 1st local day

        % Reset fields in temp
        temp.meshLightAbsorption(:) = 0;
        temp.meshThermalAbsorption(:) = 0;
        if temp.dayStep == 1
            temp.thermalDistribution(:) = 0;
        end
        
        cP = 1;  % Current partition number
        bondAlbedoArray = zeros(1, cP);

        while cP <= temp.numPartitions

            photonsPartition.values = zeros(0, 16);
            photonsPartition.thermalStart = 1;

            fprintf("\n====== PMC Partition %d of %d ======\n", ...
                    cP, temp.numPartitions)

            if size(photons.values, 1) > 0

                % Partition solar photons
                if photons.thermalStart > 1
                    % Solar photons per partition
                    sPP = PMC.solarPhotons/temp.numPartitions;
                    solarPhotonIndices = sPP*(cP-1)+1:sPP*cP;
                else
                    solarPhotonIndices = [];
                end
                solarPhotons = photons.values(solarPhotonIndices, :);
    
                % Partition thermal photons
                if temp.dayStep == 1
                    % Thermal photons per partition per facet
                    tPP = PMC.thermalPhotons/temp.numPartitions;
                    thermalPhotonIndices = zeros(1, tPP*length(sData.areaFacet));
                    for i = 1:length(sData.areaFacet)
                        start_i = photons.thermalStart + (i-1)*PMC.thermalPhotons + tPP*(cP-1);
                        end_i = start_i + tPP - 1;
                        thermalPhotonIndices((i-1)*tPP+1:i*tPP) = start_i:end_i;
                    end
                    % Every partition gets an equal number of thermal
                    % photons from each facet
                else
                    thermalPhotonIndices = [];
                end
                thermalPhotons = photons.values(thermalPhotonIndices, :);
    
                photonsPartition.values = [solarPhotons; thermalPhotons];
                photonsPartition.thermalStart = size(solarPhotons, 1) + 1;
            end  % End of partitioning prep

            % Photon tracing
            traceStart = tic();
            [aData, ~, bondAlbedo] = simulatePhotons(sData, photonsPartition, PMC);
            bondAlbedoArray(cP) = bondAlbedo;
            fprintf("\n%30s %f s\n","Photon tracing time:", toc(traceStart));

            % Binning absorption events into heat transfer mesh
            binStart = tic();
            binnedData = binAbsorptionEvents(temp, aData, energyPerFacet, cP);                
            temp.meshLightAbsorption = temp.meshLightAbsorption + binnedData.meshLightAbsorption;
            if temp.dayStep == 1
                temp.thermalDistribution = temp.thermalDistribution + binnedData.thermalDistribution;
                temp.meshThermalAbsorption = temp.meshThermalAbsorption + binnedData.meshThermalAbsorption;
            elseif cP == 1  % Prevent duplicate thermal absorption data
                temp.meshThermalAbsorption = binnedData.meshThermalAbsorption;
            end
            fprintf("\n%30s %f s\n","Binning time:", toc(binStart));

            cP = cP + 1;

            if temp.dayStep > 1 && size(photons.values, 1) == 0
                break
            end
        end  % End of partitioning while loop

        % Cache solar results of 1st local day
        temp.prevSolarAbsorption(:, temp.dayStep) = temp.meshLightAbsorption;
        temp.bondAlbedo = mean(bondAlbedoArray);
        
    else  % No partitioning after 1st local day

        % Photon tracing
        traceStart = tic();
        [aData, ~, bondAlbedo] = simulatePhotons(sData, photons, PMC);
        temp.bondAlbedo = bondAlbedo;
        fprintf("\n%30s %f s\n","Photon tracing time:", toc(traceStart));
        
        % Binning absorption events into heat transfer mesh
        binStart = tic();
        binnedData = binAbsorptionEvents(temp, aData, energyPerFacet);
        temp.meshLightAbsorption = binnedData.meshLightAbsorption;
        temp.meshThermalAbsorption = binnedData.meshThermalAbsorption;
        fprintf("\n%30s %f s\n","Binning time:", toc(binStart));

    end
            
end