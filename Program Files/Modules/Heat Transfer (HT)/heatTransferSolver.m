function temp = heatTransferSolver(temp,HT)
%       Solve for temperatures in the thermal domain based on the
%   simulation type.

fprintf("    Running heat transfer solver\n"); tic;
        % [Validation] Plot the initial temperatures
%     figure(1);
%     if temp.tInt(1) ~= 0
%         nX = temp.meshTemperatures(:,1); % [m] x coords of nodes
%         nZ = temp.meshTemperatures(:,2); % [m] z coords of nodes
%         nT = temp.meshTemperatures(:,3); % [K] old node temperatures
%         scatter(nX,nZ,60,nT,'filled');
%     end
%     colorbar; axis equal; title(sprintf("t = %f s",temp.tInt(1)));
%     drawnow;

    dt = temp.time;
    nS = length(temp.xS)-1; % number of snow surfaces

    tModel = temp.tModel;

    % Get number of timestep-independent boundary conditions
    numBC = length(tModel.BoundaryConditions.ThermalBCAssignments);

    % Assign heat fluxes to each snow surface
    if ~isempty(temp.thermalEmission)
        dA = vecnorm([diff(temp.xS);diff(temp.zS)])'; % [m^2]
        qS = temp.thermalEmission/dt./dA; % [W/m^2] fluxes
        for ii = 1:nS
            if HT.type == "Planetary"
            thermalBC(tModel,'Edge',ii+nS*[0,1,2],'HeatFlux',-qS(ii));
            elseif HT.type == "ArkChamber" || HT.type == "Open"
                thermalBC(tModel,'Edge',ii,'HeatFlux',-qS(ii));
            end
        end
    end
    
    % Assign internal heat sources
    tModel = assignHeatSources(dt,... % time step size
                               temp.meshElements,... % elements
                               temp.meshNodes,... % mesh nodes
                               temp.meshLightAbsorption,...
                               temp.meshThermalAbsorption,...
                               tModel,... % thermal model
                               HT); % Heat Transfer options

    % Assign initial temperatures
    tModel = assignTemperatures(tModel,temp.meshTemperatures,HT);
    
    % Solve for temperatures over time interval dt
    soln = solve(tModel,[0 dt]); % solve
    
    % Obtain surface temperatures
    if HT.type == "Planetary"
        sT = sampleSurfaceTemperatures(nS+(1:nS),soln);
    else
        sT = sampleSurfaceTemperatures(1:nS,soln);
    end

    % Extract solution data
    nT = soln.Temperature(:,end); % temperatures in entire mesh
    if HT.type == "Planetary"
        % Keep nodes, elements, and temperatures only for center region
        [nds,~,nId] = extractCenterMesh(tModel);
        nT = nT(nId,:); % temperatures in the center region
    else
        nds = soln.Mesh.Nodes'; % nodes
    end

    % Store data
    temp.meshTemperatures = [nds,nT];
    temp.tempProfile = sT;

    % Reset tModel by deleting timestep-dependent boundary conditions
    % Otherwise, number of BCs would increase every timestep
    tModel.BoundaryConditions.ThermalBCAssignments(numBC+1:end) = [];
    tModel.InitialConditions = [];  % Temperatures of snowpack interior
    tModel.HeatSources = [];  % Deposited solar/thermal energy

    % Print stats
    fprintf("  Avg surface temperature: %f K\n",mean(sT));
    fprintf("  Min surface temperature: %f K\n",min(sT));
    fprintf("  Max surface temperature: %f K\n",max(sT));
    fprintf("          Completion time: %f s\n\n",toc);
        % [Validation] Plot the final temperatures
%     figure(2);
%     nX = temp.meshTemperatures(:,1); % [m] x coords of nodes
%     nZ = temp.meshTemperatures(:,2); % [m] z coords of nodes
%     nT = temp.meshTemperatures(:,3); % [K] new node temperatures
%     scatter(nX,nZ,60,nT,'filled');
%     colorbar; axis equal; title(sprintf("t = %f s",temp.tInt(2)));
%     drawnow;
end