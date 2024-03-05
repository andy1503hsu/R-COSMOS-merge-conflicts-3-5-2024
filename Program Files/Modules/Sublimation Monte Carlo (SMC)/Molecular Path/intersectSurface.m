function [xSc,ind] = intersectSurface(m,n,xi,xf,zi,zf,mF,nF,zM,uC) %#codegen
%       Finds the intersection points between a line and a surface.
%
%   For the line:
%       m = slope
%       n = y-intersect
%
%   For the surface:
%       xi|xf = initial|final X coordinates
%       zi|zf = initial|final Z coordinates
%          mF = slopes
%          nF = y-intersects
%          zM = zero slopes (horizontal facets)
%          uC = undefined conditions (vertical facets)

    % Theoretical intersects
    xI = (nF-n)./(m-mF); % X coordinates
    zI = m*xI + n; % Z coordinates

    % Undefined conditions (vertical facets)
    if uC > 0
        xI(uC) = xi(uC); % X intersect coordinates
        zI(uC) = m*xI(uC) + n; % Z intersect coordinates
    end

    % Zero slopes (horizontal facets)
    if zM > 0
        zI(zM) = zi(zM); % Z coordinates (for zero-slope)
    end

    % Intersections index
    index = (xf-xI).*(xi-xI) <= 0 & (zf-zI).*(zi-zI) <= 0;
    if any(index) % for existing intersections
        ind = find(index); % index of each facet intersected by the line
        xSc = xI(ind); % X coordinates of the intersection points

    else % for no intersections
        xSc = 0; ind = 0; % return zero index
    end
end