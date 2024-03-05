function viewTemperatures(model,data,dd,ss)
%       Display the temperature distribution in the snown region over time
%   as the surface evolves.

    % Extract essential data
    ele = data.mesh.elements{dd}; % elements
    nds = data.mesh.nodes{dd}; % nodes
    nT = data.mesh.temperatures{dd}(:,end-model.resolution.daySteps+ss); % nodal temperatures
    dayN = model.resolution.maxDays; % number of days
    dayL = rotationPeriod(model.body.rotation); % [s] day length
    dayStp = model.resolution.daySteps; % time steps per day
    
    % Find domain limits and time
    xMin = min(nds(:,1)); xMax = max(nds(:,1)); % [m] domain x limits
    zMin = min(nds(:,2)); zMax = max(nds(:,2)); % [m] domain z limits
    dt = dayL/dayStp; % [s] timestep size
    tG = sub2ind([dayStp,dayN],ss,dd); % find correct time step
    tt = tG*dt; % [s] new time
    
    % Create temperature grid for contour plot
    xRes = 100; zRes = 100; % number of grid points in x and z direction
    [xG,zG] = ndgrid(linspace(xMin,xMax,xRes),linspace(zMin,zMax,zRes));
    TG = griddata(nds(:,1),nds(:,2),nT,xG,zG);
    msh = triangulation(ele,nds);
    TG(isnan(pointLocation(msh,xG(:),zG(:)))) = nan;
    
    % Create initial plots
    P = patch('Faces',ele,'Vertices',nds,'FaceVertexCData',nT); hold on;
    %contour(xG,zG,TG,'LineColor','k');
    B = colorbar; hold off;
    
    % Adjust and label
    P.FaceColor = 'interp'; P.EdgeColor = 'none';
    xlim([xMin xMax]); ylim([zMin zMax]);
    title(sprintf("t = %0.2f s",tt));
    B.Label.String = "Temperature [K]";
    xlabel("X [m]"); ylabel("Z [m]");
    colormap('Jet');
end