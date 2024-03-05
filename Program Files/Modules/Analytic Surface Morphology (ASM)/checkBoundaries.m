function [endX,endY] = checkBoundaries(Surface,finalX,finalY) %#codegen
%       Ensures that the boundary conditions are met according to the
%   simulation type.
    
    % Extract the final surface information (only interested in the facets)
    TempFacets = extractSurfaceData(finalX,finalY);
    
    % Create an index to mark the location of surface crossings
    idx = false(1,length(TempFacets.xi));
    
    % Determine crossings between facets near the boundaries
    idx = locateSurfaceCrossings(finalX,finalY,idx);
    
    % Remove the marked facets
    TempFacets = removeMarkedFacets(TempFacets,idx);
    
    % Crop intersecting facets to coincide on the intersection point
    TempFacets = intersectFacets(TempFacets);
    
    % Extract the final nodes from the temporary facets structure
    endX = [TempFacets.xi TempFacets.xf(end)];
    endY = [TempFacets.yi TempFacets.yf(end)];
    
    % Boundary conditions test
    if ~any(endX < Surface.xS(1) | endX > Surface.xS(end))
        fprintf('  Boundary conditions test: passed\n');
    else
        warning('Check boundary conditions')
    end
end