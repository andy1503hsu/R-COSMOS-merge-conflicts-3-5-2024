function viewLightSourcePath(model,oS)


    % Generate orbital path
    orbit = orbitalPath(model);

    % Generate light source path
    light = lightPath(orbit,model);
% light =path
    plot3(light.lightPath(1,:,oS),light.lightPath(2,:,oS),light.lightPath(3,:,oS),'bo')
    hold on

    % Generate surface
    surface = surfaceShape(model.surface);
    surface.z = (surface.z - surface.depth)/max(surface.z) - 0.5;
    surface.x = surface.x/max(surface.x) - 0.5;
    
    plot3(surface.x,zeros(1,length(surface.x)),surface.z)
    
    axis equal
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
    xlabel('x')
    ylabel('y')
    zlabel('z')

end

