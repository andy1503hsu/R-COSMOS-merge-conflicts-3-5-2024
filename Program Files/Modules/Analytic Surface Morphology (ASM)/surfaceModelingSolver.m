function Surface = surfaceModelingSolver(Surface,displacement,Opt) %#codegen
%       Given a surface shape and the magnitude of its displacement, this
%   function reconstructs the surface and returns a new surface.
    
    Surface.yS = Surface.zS; % temporary correction for coordinate variable

    fprintf('Analytic Surface Morphology (Ver. 2.1)\n'); startT = tic;
    fprintf('      Average displacement: %.2e m\n',mean(displacement));
    
    % 1. Create the simulation model according to the user-selected options
    Opt = checkUserInputs(Opt);
    
    % 2. Estimate the precision of the reconstruction
    modelingPrecision(Surface,displacement);
    
    % 3. Discretize the surface for the selected simulation type
    [xS,yS,dS,ind] = prepareSurface(Surface,displacement,Opt);
    
    % 4. Extract the surface information (area, slopes, y-intersects, ...)
    sData = extractSurfaceData(xS,yS);
    
    % 5. Find the location of critical points in the original surface
    [fLoc,nLoc] = findCriticalLocations(sData,dS,Opt);
    
    % 6. Calculate the position of each facet in space and the travel path
    % of adjacent surface nodes
    [Facets,Nodes] = displaceFacets(sData,dS);
    
    % 7. Resize the facets to reach the path of adjacent surface nodes
%     Facets = resizeFacets(Facets,Nodes,nLoc,dS);
    
    % 8. Adjust the surface at the critical locations (if desired)
    if Opt.declare.adjustByDistance
        ind = adjustCriticalFacets(fLoc,sData,Nodes,dS,ind);
    end
    
    % 9. Average the endpoints of adjacent facets to create temporary nodes
    [tempX,tempY] = reconstructSurface(Facets);
    
    % 10. Find and mark the locations of facets in between crossings
    ind = locateSurfaceCrossings(tempX,tempY,ind);
    
    % 11. Remove the marked facets
    Facets = removeMarkedFacets(Facets,ind);
    
    % 12. Crop intersecting facets to coincide on the intersection point
%     Facets = intersectFacets(Facets);
    
    % 13. Average the endpoints of adjacent facets to create new nodes
    [midX,midY] = reconstructSurface(Facets);
    
    % 14. Revise the surface for intersections created by step 12
    [Facets,midX,midY] = reviseSurface(Facets,midX,midY,Opt);
    
    % 15. Apply boundary conditions
    [finalX,finalY] = boundaryConditions(Surface,Facets,midX,midY,Opt);
    
    % 16. Check the boundary conditions for irregularities
    [endX,endY] = checkBoundaries(Surface,finalX,finalY);
    
    % 17. Restore compatibility
    Surface = restoreCompatibility(Surface,endX,endY);
    
    % 18. Determine surface stability
    testSurfaceStability(Surface)
    
    fprintf('           Simulation time: %.2f seconds\n',toc(startT));

    Surface.zS = Surface.yS; % temporary correction for coordinate variable
end