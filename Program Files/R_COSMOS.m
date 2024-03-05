function data = R_COSMOS(model)
%       Simulates the surface morphology over time of a given material in
%   airless environments subject to light and thermal radiation, conductive
%   heat transfer, and sublimation. The results are simulated given an
%   initial surface shape and material properties.

fprintf('R-COSMOS ver. 5.3\n')

% Setup the R-COSMOS simulation environment
[model, temp, data, PMC, HT, SMC, ASM] = modelSetup(model);

% Start the simulation
startT = tic; % simulation time
for gS = 1:model.resolution.geologicSteps % for each geologic period
    temp.geoStep = gS; % update the geologic step

    for oS = 1:model.resolution.orbitSteps % for each orbital period
        temp.orbitStep = oS; % update the orbital step
        temp.lightDistance = data.orbit.radii(oS); % [AU] light distance
        
        % Surface discretization (model architecture)
        sData = surfaceData(temp);

        % Simulate the diurnal cycle until it reaches steady-state
        temp.iter = 1;
        while true % diurnal cycle
            % Restart head
            if temp.iter > 1
                % Clear the cashed information
                clear all
    
                % Restart the simulation
                load Temporary/restartHead/restartData
            end
    
            % Radiative Heat Transfer process
            for dS = 1:model.resolution.daySteps % for each timestep
                temp.dayStep = dS; % update the day step

                % Update zenith and azimuth
                temp.zen = data.light.zenith(dS:dS+1,oS); % [deg]
                temp.azi = data.light.azimuth(dS:dS+1,oS); % [deg]
                temp.lightPath = data.light.lightPath(:,dS:dS+1,oS);
                
                % Display the status of the simulation
                simulationStatus(temp, model, 'RHT')
 
                % Photon Monte Carlo (solar and thermal radiation)
                temp = radiationSolver(temp, sData, PMC);
    
                % Heat transfer
                temp = heatTransferSolver(temp, HT);
    
                % Data storage
                data = dataStorage(temp, data, 'RHT');
            end
    
            %parfor iT = 1:3
            %  disp(iT)
            %end
            % Sublimation process
            % If the diurnal cycle has reached steady-state
            if stoppingCondition(temp, data, model, 'steadyState')
                for dS = 1:model.resolution.daySteps % for each timestep
                    temp.dayStep = dS; % update the local step
    
                    % Update zenith and azimuth, needed for solar exitance
                    temp.zen = data.light.zenith(dS:dS+1,oS); % [deg]
                    temp.azi = data.light.azimuth(dS:dS+1,oS); % [deg]

                    % Display the status of the simulation
                    simulationStatus(temp, model, 'SMC')
    
                    % Recall surface temperature profile
                    temp.tempProfile = data.surface.temperature(dS,:,oS,gS);
               
                    % Get thermal exitance results
                    temp = getThermalExitance(temp, sData, PMC);

                    % Get solar exitance results
                    temp = getSolarExitance(temp, sData, PMC);

                    % Sublimation Monte Carlo (sublimation)
                    temp = sublimationSolver(temp,SMC);
    
                    % Data storage
                    data = dataStorage(temp, data, 'SMC');
                end
                break % since the diurnal cycle is at steady-state
    
            else % if not at steady-state
                temp.iter = temp.iter + 1; % update the iteration count
            end
    
            % Update the surface temperature
            temp.actualSurfTemp = data.surface.temperature(:,:,oS,gS);
    
            % Save the simulation progress
            save('Temporary/restartHead/restartData', "-v7.3")
        end

        % Surface morphology changes along the orbital path
        if model.orbitalMorphology
            % Net displacement of each surface facet
            temp = projectedDisplacement(temp, model, data, 'orbital');

            % Surface morphology
            temp = surfaceModelingSolver(temp, temp.netDisplacement, ASM);

            % Update thermal model
            temp = thermalModel(temp,HT);

            % Data storage
            data = dataStorage(temp, data, ['ASM', 'orbital']);
        end
    end

    % Surface morphology changes at a geologic scale
    if ~model.orbitalMorphology
        % Net displacement of each surface facet
        temp = projectedDisplacement(temp, model, data, 'geological');
    
        % Surface morphology
        temp = surfaceModelingSolver(temp, temp.netDisplacement, ASM);

        
        %[temp.xS,temp.zS] = surfaceModelling(temp,temp.netDisplacement,model.type);
        
        % Update thermal model
        temp = thermalModel(temp,HT);

        % Data storage
        data = dataStorage(temp, data, ['ASM', 'geological']);
    end

    % Save the simulation progress
    save('Cases/' + model.filename + ".mat", "model", "data", "-v7.3")
end

% Total simulation time
fprintf("\n Total completion time: "); timeDisplay(toc(startT))

% Export R-COSMOS log
diary off
end