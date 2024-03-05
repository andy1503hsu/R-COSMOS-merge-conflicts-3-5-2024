function [xi,xf,zi,zf,mF,nF,zM,uC] = commonVariables(sData) %#codegen
%       Selects the variables to be used by the ray-trace algorithm and
%   identifies the facets with zero slopes and/or undefined conditions.

    % Surface information
    xi = sData.xi; xf = sData.xf; % initial|final X coordinates
    zi = sData.zi; zf = sData.zf; % initial|final Y coordinates
    
    mF = sData.mFacet; % slopes of the surface facets
    nF = sData.nFacet; % y-intersects of the surface facets

    % Zero slopes (horizontal facets)
    if any(mF == 0) % for any facet with slope equal to zero
        zM = find(mF == 0); % index of the zero-slope facets
    else
        zM = 0; % variable set to 0
    end
    
    % Undefined conditions (vertical facets)
    if any(mF == Inf)% for any facet with undefined conditions
        uC = find(mF == Inf); % index of the undefined facets
    else
        uC = 0; % variable set to 0
    end
end