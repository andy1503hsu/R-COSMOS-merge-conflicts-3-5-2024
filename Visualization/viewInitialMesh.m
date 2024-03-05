function viewInitialMesh(model)
%   Visualize the starting mesh.
    
    % Extract essential data
    xS = model.surface.x;
    zS = model.surface.z;
    Hmin = model.resolution.Hmin;
    Hmax = model.resolution.Hmax;
    Hgrad = model.resolution.Hgrad;

    if model.type == "Planetary"
        msh = meshPlanetary(xS,zS,Hmin,Hmax,Hgrad);
        
    elseif model.type == "ArkChamber"
        msh = meshArkChamber(xS,zS,Hmin,Hmax,Hgrad);
        
    end
    
    % Make initial plots
    triplot(msh); hold on;
    plot(xS,zS,'k','LineWidth',2);
    scatter(xS,zS,30,'k','filled');
    
    % Label and adjust
    title(sprintf("Hmin = %0.4f, Hmax = %0.4f, Hgrad = %0.4f",...
          Hmin,Hmax,Hgrad));
    xlabel("X [m]"); ylabel("Y [m]");
    axis equal
    
end