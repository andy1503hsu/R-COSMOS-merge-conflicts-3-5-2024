function [xSc,zSc] = intersectSurfacePMC(x1,z1,x2,z2,xi,zi,xf,zf,...
                                         zM,uC,x34,z34,determinant34) %#codegen
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
    [pX,pZ] = intersectLinesPMC(x1,z1,x2,z2,x34,z34,determinant34);

    % Take care of edge cases (horizontal and vertical facets)
    
    if zM > 0 % zero-slope facets exist
        pZ(zM) = zi(zM); % z-coord of intersection will just be z-coord of that facet
    end
    
    if uC > 0 % vertical facets exist
        pX(uC) = xi(uC); % x-coord of intersection will just be x-coord of that facet
        
        % Line properties
        m = (z1-z2)/(x1-x2); % slopes
        n = z1 - m.*x1; % z-intersects
        pZ(uC) = m*pX(uC) + n; % Z intersect coordinates
    end
    
    % Intersections index
    index = (xf-pX).*(xi-pX) <= 0 & (zf-pZ).*(zi-pZ) <= 0;
    if any(index) % for existing intersections
        ind = find(index); % index of each facet intersected by the line
        xSc = pX(ind); % X coordinates of the intersection points
        zSc = pZ(ind);
    else % for no intersections
        xSc = 0; zSc = 0;  % return zero index
    end

    %{
    % Slower 
    xSc = pX(index);
    if isempty(xSc)
        xSc = 0; zSc = 0;
    else
        zSc = pZ(index);
    end
    %}
end