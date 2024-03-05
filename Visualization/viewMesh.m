function viewMesh(model,data)
%       Display the thermal finite element mesh of the snow region over
%   time as the surface evolves.

    % Extract initial data to plot
    ele = data.mesh.elements{1}; % elements
    nds = data.mesh.nodes{1}; % nodes
    xS = model.surface.x; % [m] x coords of surface nodes
    zS = model.surface.z; % [m] z coords of surface nodes
    wD = max([max(data.surface.x,[],'all'),xS]); % [m] domain width
    hD = max([max(data.surface.z,[],'all'),zS]); % [m] domain height
    
    % Create initial plots
    figure();
    P = patch('Faces',ele,'Vertices',nds); hold on;
    N = plot(xS,zS,'k','LineWidth',2);
    S = scatter(xS,zS,20,'k','Filled');
    
    % Adjust and label
    P.FaceColor = 'none';
    %axis equal; 
    xlim([0 wD]); ylim([0 hD]);
    title(sprintf("t = %.2f s (%.2f days)",0,0));
    xlabel("X [m]"); ylabel("Z [m]");
    
    for day = 1:model.time.days
        
        
        % Extract new surface and mesh
        if day == 1
            xS = model.surface.x; % [m] x coords of surface nodes
            zS = model.surface.z; % [m] z coords of surface nodes
            
        else
            xS = data.surface.x(day-1,:); % [m] x coords of surface nodes
            zS = data.surface.z(day-1,:); % [m] z coords of surface nodes
        end
        
        ele = data.mesh.elements{day}; % elements
        nds = data.mesh.nodes{day}; % nodes
        
        % Update graphics
        P.Faces = ele;
        P.Vertices = nds;
        N.XData = xS;
        N.YData = zS;
        S.XData = xS;
        S.YData = zS;

        title(sprintf("day %d",day));
        
        if day ~= model.time.days
            pause()
        end
    end
    
end