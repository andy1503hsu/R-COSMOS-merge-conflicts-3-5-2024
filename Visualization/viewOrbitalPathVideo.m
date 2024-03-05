function viewOrbitalPathVideo(model)


    % True anomaly
    iTa = model.orbit.initialTrueAnomaly; % [deg]
    theta = wrapTo360(linspace(iTa, 360 + iTa, model.resolution.orbitSteps+1)); % [deg]

    % Generate orbital path
    orbit = orbitalPath(model);

    % Generate light source path
    light = lightPath(orbit,model);

    % Generate surface
    surface = surfaceShape(model.surface);

    % Scale the surface to unity and center
    surface.z = (surface.z - surface.depth)/max(surface.z) - 0.5;
    surface.x = surface.x/max(surface.x) - 0.5;

    % Orbit
    subplot(1,2,1)
    polarplot(deg2rad(theta),orbit.radii,...
            'LineWidth',1.5,...
            'Color',"#A2142F")
    hold on
    polarplot(0,0,...
            'Marker','o',...
            'MarkerSize',40,...
            'MarkerFaceColor',"#EDB120",...
            'MarkerEdgeColor','black')
    title('Orbital path (in AU)')

    % Light path
    subplot(1,2,2)
    plot3(surface.x,zeros(1,length(surface.x)),surface.z,...
        'black',...
        'LineWidth',1.5)
    hold on
    title('Light Source Path')
    axis equal
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
    xlabel('x (m)','FontSize',12)
    ylabel('y (m)','FontSize',12)
    zlabel('z (m)','FontSize',12)

    % Coordinate axis
    % quiver3(0,0,-0.5,0.25,0,0,'k','LineWidth',1,'AutoScale','on')
    % quiver3(0,0,-0.5,0,0.25,0,'k','LineWidth',1,'AutoScale','on')
    % quiver3(0,0,-0.5,0,0,0.25,'k','LineWidth',1,'AutoScale','on')
for kk = 1:2
    for ii = 1:length(theta)
        subplot(1,2,1)
        h = polarplot(theta(ii)*pi/180,orbit.radii(ii),...
                    'Marker','o',...
                    'MarkerSize',20,...
                    'MarkerFaceColor',"#0072BD",...
                    'MarkerEdgeColor','black');

        subplot(1,2,2)
        g = plot3(light.lightPath(1,:,ii),light.lightPath(2,:,ii),light.lightPath(3,:,ii),...
                'Color',"#0072BD",...
                'LineWidth',1.5);
        for jj = 1:length(light.lightPath(1,:,1))-1
            f = scatter3(light.lightPath(1,jj,ii),light.lightPath(2,jj,ii),light.lightPath(3,jj,ii),...
                100,...
                [0.9290 0.6940 0.1250],...
                'filled');
            pause(0.01)
            delete(f)
        end
% % pause(0.1)
        delete(h)
        delete(g)
    end
end
    % Final plots
    % Orbit
    subplot(1,2,1)
    polarplot(deg2rad(theta),orbit.radii,...
            'LineWidth',1.5,...
            'Color',"#A2142F")
    % Light path
    subplot(1,2,2)
    plot3(light.lightPath(1,:,ii),light.lightPath(2,:,ii),light.lightPath(3,:,ii),...
            'Color',"#0072BD",...
            'LineWidth',1.5);

end