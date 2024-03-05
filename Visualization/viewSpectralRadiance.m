function viewSpectralRadiance(model,data,dd,ss,wId,sId,src)
%       View the spectral radiance at the specified day and time at the
%   specified surface, in the specified wavelength interval.
    
    % Extract the correct exitance energy grid
    if lower(src) == "light"
        exiG = squeeze(data.energy.light.exitance(dd,ss,wId,sId,:,:)); % [J]
    elseif lower(src) == "thermal"
        exiG = squeeze(data.energy.thermal.exitance(dd,ss,wId,sId,:,:)); % [J]
    else
        error("Unknown light source type");
    end
    
    % Create zenith and azimuth bin edges
    zenN = model.resolution.zenithBins; % number of zenith bins
    aziN = model.resolution.azimuthBins; % number of azimuth bins
    zenE = linspace(0,pi/2,zenN+1); % [rad] zenith bin edges
    aziE = linspace(0,2*pi,aziN+1); % [rad] azimuth bin edges
    
    % Dimension interval sizes for each dimension
    dt = model.time.dayLength/model.time.daySteps; % [s] time step size
    dw = diff(model.light.wavelengthBounds)/...
         model.resolution.wavelengthBins; % [m] wavelength interval size
    dA = max(model.surface.x)-min(model.surface.x)/...
         model.resolution.ceilingSegments; % [m^2] exitance area
    
    % Obtain the spectral irradiance and other related quantities
    [sR,~,~,pX,pY] = getSpectralRadiance(dt,dw,dA,zenE,aziE,exiG);
    
    % Find the correct time
    tG = sub2ind([model.time.daySteps,model.time.days],ss,dd);
    tt = tG*dt; % [s] time
    
    % Find the correct wavelength interval
    wE = linspace(model.light.wavelengthBounds(1),...
                  model.light.wavelengthBounds(2),...
                  model.resolution.wavelengthBins+1);
    wL = wE(wId); wU = wE(wId+1);
    
    % Create initial plots
    figure();
    sR_ = zeros(size(sR)+1);
    sR_(1:end-1,1:end-1) = sR;
    P = pcolor(pX,pY,sR_);
    C = colorbar;
    
    % Adjust and label
    P.EdgeColor = "None";
    colormap('jet');
    view(0,90); axis equal;
    title(sprintf('t \\in [%f %f] s, \\lambda \\in [%f %f] \\mum',...
                        tt-dt,tt,wL*1e6,wU*1e6));
    C.Label.String = "Spectral Radiance [W/m^2/sr/m]";
    
end