function viewSurfaceEvolution(model,data)
%       Visualizes the surface evolution over time.

    % Initial surface
    Color = [0 71/255 100/255]; % neutral surface color
    plot(model.surface.x,model.surface.z,'-','LineWidth',3,'Color',Color)
% pause(0.5)
    hold on
    % Surface evolution over time
    for day = 1:model.resolution.geologicSteps
        plot(data.surface.x(1,:,day),data.surface.z(1,:,day),'b-')
%         plot(data.surface.x(day,:),data.surface.z(day,:),'ro')
        axis equal
        hold on
        if data.surface.z(1,:,day) == 0
            break
        end
%         pause(0.1)
    end
    
%     figure()
%     for day = 1:model.time.days
%         plot(day,max(data.surface.z(day,:))-min(data.surface.z(day,:)),'bo');
%         hold on
%     end
end