function viewZenithAzimuthPath(model,tL,tU)
%       Display the path of the light source for planetary type simulation
%   in a celestial sphere in the time interval specified. The domain is 
%   centered in this celestial sphere.

    if lower(model.type) == "ark chamber"
        error("Simulation type must be planetary");
    end
    
    % Extract the zenith and azimuth path functions
    zenF = model.light.zenithPath;
    aziF = model.light.azimuthPath;
    
    % Create a parametric line of the light source path
    rxL = @(t) sin(zenF(t)).*cos(aziF(t));
    ryL = @(t) sin(zenF(t)).*sin(aziF(t));
    rzL = @(t) cos(zenF(t));
    
    % Create a parametric line of the celestial horizon
    rxH = @(t) cos(t);
    ryH = @(t) sin(t);
    rzH = @(t) zeros(size(t));
    
    % Create initial plots
    fplot3(rxL,ryL,rzL,[tL,model.time.dayLength],'--','LineWidth',2,'Color','#A2142F'); hold on;
%     fplot3(rxH,ryH,rzH,[0 2*pi],'k','LineWidth',2);
%     scatter3(rxL(tL),ryL(tL),rzL(tL),80,'filled','MarkerEdgeColor','#A2142F','MarkerFaceColor', '#A2142F');
    scatter3(rxL(tU),ryL(tU),rzL(tU),160,'filled','MarkerEdgeColor', '#A2142F', 'MarkerFaceColor','white','LineWidth',2);
    
    % Adjust and Label
    A = gca;
    axis equal;
    xlim([-1 1]); ylim([-1 1]); zlim([-1 1]); 
    xlabel("X"); ylabel("Y"); zlabel("Z");
    A.XAxisLocation = 'origin';
    A.YAxisLocation = 'origin';
%     set(gca,'visibility','off');
end