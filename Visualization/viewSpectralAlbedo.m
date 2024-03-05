function [wE,alb] = viewSpectralAlbedo(model,data,dd,ss)
%       Display the temperature distribution in the snown region over time
%   as the surface evolves.

    % Extract and compute essentials
    wBins = model.resolution.wavelengthBins; % number of wavelength bins
    wLim = model.light.wavelengthBounds; % [m] wav. int.
    emiG = data.energy.light.emission; % [J] emission grid
    exiG = data.energy.light.exitance; % [J] exitance grid
    
    dt = model.time.dayLength/model.time.daySteps; % [s] time step size
    wE = linspace(wLim(1),wLim(2),wBins+1); % [m] wavelength bin edges
    
    tG = sub2ind([model.time.daySteps,model.time.days],ss,dd);
    tt = tG*dt; % [s] update time

    % Compute spectral albedo
    eIn = emiG(dd,ss,:); % [J] incident energy
    eOut = sum(exiG(dd,ss,:,:,:,:),[4 5 6]); % [J] reflected energy
    alb = squeeze(eOut)./squeeze(eIn);
    alb(isnan(alb)) = 0;
    
    % Create initial plot
    H = stairs(wE,[alb;alb(end)]);
    
    % Adjust and label
    H.LineWidth = 2; H.Color = 'k';
    
    ylim([0 1]);
    xlabel("Wavelength [m]");
    ylabel("Albedo [-]");
    title(sprintf("t \\in [%.2f %.2f] s",tt-dt,tt));
    set(gca,'FontName','Calibri','FontSize',12);
end