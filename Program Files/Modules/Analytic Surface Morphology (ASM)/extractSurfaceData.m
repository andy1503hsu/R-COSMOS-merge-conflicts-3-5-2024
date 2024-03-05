function sData = extractSurfaceData(xS,yS) %#codegen
%       Extracts the surface information based on the given shape
%   distribution and stores it in a structure. Establishes the definitions
%   of ASM for undefined conditions.

    % Extract the surface facets coordinates
    [xi,xf,yi,yf] = discretizeSurface(xS,yS);
    
    %    Initial                      Final
    sData.xi = xi;           sData.xf = xf; % X coordinates [m]
    sData.yi = yi;           sData.yf = yf; % Y coordinates [m]

    % Area of the facets (length for 2D model)
    areaFacet = calculateFacetsArea(xi,xf,yi,yf);
    sData.areaFacet = areaFacet; % [m^2]

    % Slopes and y-intercepts of the facets (formula for a line)
    %    Formulas                     ASM definitions
    mF = (yf-yi)./(xf-xi);      mF(mF == -Inf) = Inf; % slopes
    nF = yf - mF.*xf;           nF(isnan(nF)) = -Inf; % y-intercepts

    % Storage
    sData.mFacet = mF; % final slopes
    sData.nFacet = nF; % final y-intercepts
    
    % Unit normal vector
    angleFacet = acosd((xf-xi)./areaFacet); % [deg]
    % Check and adjust the angle to point outwards from the surface
    for ind = 1:length(xi)
        if yf(ind) > yi(ind)
            angleFacet(ind) = 90 + angleFacet(ind); % [deg]
        elseif yf(ind) < yi(ind)
            angleFacet(ind) = 90 - angleFacet(ind); % [deg]
        else
            angleFacet(ind) = 90 + 0; % [deg]
        end
    end
    sData.angleFacet = angleFacet; % [deg]
end