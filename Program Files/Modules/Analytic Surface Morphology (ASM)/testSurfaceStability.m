function testSurfaceStability(Surface) %#codegen
%       Determine if the final surfaces is stable, hence free of
%   intersections between its facets.

    % Extract the final surface information (only interested in the facets)
    TempFacets = extractSurfaceData(Surface.xS,Surface.yS);
    
    % Create an index to store the location of possible intersections
    idx = false(1,length(TempFacets.xi));
    
    % Check for any intersections between the facets
    idx = locateSurfaceCrossings(Surface.xS,Surface.yS,idx);
    
    % Print the results of the stability test
    if any(idx) % if there is any intersection the surface is not stable
        warning('Surface is not stable.')
    else % the surface is stable
        fprintf('            Stability test: passed\n');
    end
end