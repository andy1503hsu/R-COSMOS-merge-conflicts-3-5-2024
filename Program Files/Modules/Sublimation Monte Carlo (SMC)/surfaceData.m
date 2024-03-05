function sData = surfaceData(in) %#codegen
%       Extracts the surface information based on the given shape
%   distribution and stores it in a structure. Establishes the definitions
%   of SMC for undefined conditions.

    % Surface coordinates
    %    Initial                      Final
    xi = in.xS(1:end-1);        xf = in.xS(2:end); % temporary X [m]
    zi = in.zS(1:end-1);        zf = in.zS(2:end); % temporary Z [m]
    sData.xi = xi;              sData.xf = xf; % X coordinates [m]
    sData.zi = zi;              sData.zf = zf; % Z coordinates [m]

    % Area of the facets (length for 2D model)
    areaFacet = ((xf-xi).^2 + (zf-zi).^2).^0.5; % distance formula
    sData.areaFacet = areaFacet; % [m^2]

    % Slopes and y-intercepts of the facets (formula for a line)
    %    Formulas                     SMC definitions
    mF = (zf-zi)./(xf-xi);      mF(mF == -Inf) = Inf; % slopes
    nF = zf - mF.*xf;           nF(isnan(nF)) = -Inf; % y-intercepts

    % Storage
    sData.mFacet = mF; % final slopes
    sData.nFacet = nF; % final y-intercepts
    
    % Unit normal vector
    angleFacet = acosd((xf-xi)./areaFacet); % [deg]
    % Check and adjust the angle to point outwards from the surface
    for ind = 1:length(xi)
        if zf(ind) > zi(ind)
            angleFacet(ind) = 90 + angleFacet(ind); % [deg]
        elseif zf(ind) < zi(ind)
            angleFacet(ind) = 90 - angleFacet(ind); % [deg]
        else
            angleFacet(ind) = 90 + 0; % [deg]
        end
    end
    sData.angleFacet = angleFacet; % [deg]
end