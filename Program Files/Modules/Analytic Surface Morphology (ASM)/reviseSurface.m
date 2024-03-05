function [Facets,mX,mY] = reviseSurface(Facets,mX,mY,Opt) %#codegen
%       Revises the surface for crossings of surface facets created by the
%   resizing process.
    
    % Create an index to mark the facets for elimination
    idx = false(1,length(Facets.xi));
    
    % Determine if there are crossings between facets
    idx = locateSurfaceCrossings(mX,mY,idx);
    
    % Run until the surface is stable (not common, usually only one
    % iteration is required)
    iter = 0; % start the iteration count
    while any(idx) || iter > Opt.maxiter
        
        % Remove the marked facets
        Facets = removeMarkedFacets(Facets,idx);
        
        % Crop intersecting facets to coincide on the intersection point
        Facets = intersectFacets(Facets);
        
        % Average the endpoints of adjacent facets
        [mX,mY] = reconstructSurface(Facets);
        
        % Recreate the index array to test for crossings to ensure that the
        % mid surface is stable
        idx = false(1,length(Facets.xi));
        
        % Determine if there are crossings between facets
        idx = locateSurfaceCrossings(mX,mY,idx);
        
        % Number of iterations completed
        iter = iter + 1;
    end
end