function viewIrradiance(model,data,src,day,stp,wav,bnd,zen,azi)
%       Display an image of the ceiling or floor boundary from the
%   viewing zenith and azimuth angle specified in the wavelength band
%   specified at the specified time and from the specified light source

    % Extract initial data
    xE = model.surface.x; % x edges of the ceiling boundary segments
    dA = xE(2:end) - xE(1:end-1); % [m^2] area of boundary segments
    dt = model.time.dayLength/model.time.daySteps; % [s] time step size
    
    % Grab the appropriate energy data
    if lower(bnd) == "ceiling" && lower(src) == "light"
        val = data.energy.light.exitance(day,stp,wav,2:end,zen,azi);
    elseif lower(bnd) == "floor" && lower(src) == "light"
        val = data.energy.light.exitance(day,stp,wav,1,zen,azi);
    elseif lower(bnd) == "ceiling" && lower(src) == "thermal"
        val = data.energy.thermal.exitance(day,stp,wav,2:end,zen,azi);
    elseif lower(bnd) == "floor" && lower(src) == "thermal"
        val = data.energy.thermal.exitance(day,stp,wav,1,zen,azi);
    else
        error("Unknown boundary and light source specified.");
    end
    
    % Prepare energy data to be displayed as irradiance
    val = squeeze(sum(val,[3,5,6]))'/dt./dA; % [W/m^2]
    
    % Min and max x and y values for scaling
    yMin = 0; yMax = max(val);
    xMin = min(xE); xMax = max(xE);
    
    % Plot data
    figure(); hold on;
    H = histogram('BinEdges',xE,'BinCounts',val,'DisplayStyle','stairs');
    
    % Adjust and label
    H.LineWidth = 2; H.EdgeColor = 'k';
    xlim([xMin xMax]); ylim([yMin yMax]);
    xlabel("Position [m]"); ylabel("Irradiance [W/m^2]");
    
end