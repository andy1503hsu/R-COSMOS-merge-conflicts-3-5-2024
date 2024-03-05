function viewSurfaceTemperature(model,data)
%This function plots the surface temperatures of a given model
%   outpus a figure of surface temperature
% figure ()
for ii = 1:model.resolution.geologicSteps          % number of days
    for jj = 1:model.resolution.daySteps     % number of steps
        plot(data.surface.temperature(jj,: ,ii)); %hold on;
        ylim([80 140])
        pause(0.6)
        
    end
end

end

