function data = dataStorage(temp, data, marker) %#codegen
%       Store data from each module permanently from the temporary output
%   data structure

    gS = temp.geoStep;  % Current global step
    oS = temp.orbitStep;   % Current orbit step
    iter = temp.iter;  % Current local step
    dS = temp.dayStep;  % Current step of day

    if strcmp(marker,'SMC')
        % Surface data
        data.surface.dispFacet(dS,:,oS,gS) = temp.dispFacet;
        data.surface.netFlux(dS,:,oS,gS) = temp.netFlux;

        % Exitance data
        data.exitance.solar{dS, oS, gS} = temp.solarExitance;
        data.exitance.thermal{dS, oS, gS} = temp.thermalExitance;


    elseif strcmp(marker,'RHT')
        
        if iter == 1

           if dS == 1  % Initialize mesh-dependent variables
                % Mesh data
                data.mesh.nodes{oS, gS} = temp.meshNodes;
                data.mesh.elements{oS, gS} = temp.meshElements;

                % Mesh properties
                numNodes = length(temp.meshNodes);
                numElements = length(temp.meshElements);
                numDaySteps = size(temp.prevSolarAbsorption, 2);

                % Pre-allocate mesh-dependent arrays
                data.mesh.temperatures{oS, gS} = zeros(numNodes, numDaySteps);
                data.mesh.lightAbsorption{oS, gS} = zeros(numElements, numDaySteps);
                data.mesh.thermalAbsorption{oS, gS} = zeros(numElements, numDaySteps);
           end

           % Solar photon tracing is only run on the 1st day
           data.mesh.lightAbsorption{oS, gS}(:, dS) = temp.meshLightAbsorption;

        end
    
        % Only save steady-state temperatures and thermal absorption
        % Data will be continuously overwritten until the steady-state day
        data.mesh.temperatures{oS, gS}(:, dS) = temp.meshTemperatures(:, 3);
        data.mesh.thermalAbsorption{oS, gS}(:, dS) = temp.meshThermalAbsorption;

        % From Antonio:
        data.surface.temperature(dS,:,oS,gS) = temp.tempProfile;
        data.surface.steadyStateT(dS,:,iter,oS,gS) = temp.tempProfile;

    elseif strcmp(marker,['ASM', 'orbital'])
        % Store the new surface shape
        data.surface.x(oS,:,gS) = temp.xS; % x coordinate
        data.surface.z(oS,:,gS) = temp.zS; % z coordinate

    elseif strcmp(marker,['ASM', 'geological'])
        % Store the new surface shape
        for ii = 1:oS
            data.surface.x(oS,:,gS) = temp.xS; % x coordinate
            data.surface.z(oS,:,gS) = temp.zS; % z coordinate
        end
    else
        error('Marker not found: select 1) SMC, 2) RHT, or 3) ASM')

    end
end