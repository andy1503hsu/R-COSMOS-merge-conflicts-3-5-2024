function viewHeatSources(model,data,source,dd,ss)
%       Display the energy absorbed per unit volume per unit time in the
%   mesh elements over time.
    
    % Grab the energy absorption from the specified radiation source
    if lower(source) == "light"
        eA = data.mesh.lightAbsorption{dd}(:,ss);
    elseif lower(source) == "thermal"
        eA = data.mesh.thermalAbsorption{dd}(:,ss);
    else
        error("Unkown radiation source specified");
    end
    
    % Extract and compute essential data
    ele = data.mesh.elements{dd}; % elements
    nds = data.mesh.nodes{dd}; % nodes
    dt = model.time.dayLength/model.time.daySteps; % [s] time step size
    eQ = getHeatSources(dt,ele,nds,eA);
    tG = sub2ind([model.time.daySteps,model.time.days],ss,dd);
    tt = tG*dt; % [s] time
    xMin = min(nds(:,1)); xMax = max(nds(:,1)); % [m] domain x limits
    zMin = min(nds(:,2)); zMax = max(nds(:,2)); % [m] domain z limits
    
    % Create initial plots
    hold on;
    P = patch('Faces',ele,'Vertices',nds,'FaceVertexCData',eQ);
    C = colorbar;
    
    % Adjust and label
    xlim([xMin xMax]); ylim([zMin zMax]); axis equal;
    P.FaceColor = 'flat'; P.EdgeColor = 'none';
    C.Label.String = "Heat Source [W/m^3]";
    xlabel("X [m]"); ylabel("Z [m]");
    title(sprintf("t \\in [%.2f %.2f] s",tt-dt,tt));
end