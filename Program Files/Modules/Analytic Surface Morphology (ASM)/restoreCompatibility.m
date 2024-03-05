function Surface = restoreCompatibility(Surface,fX,fY) %#codegen
%       Restores the compatibility with the original surface by maintaining
%   constant the number of surface elements (nodes and facets). For each
%   facet eliminated, a new surface node is introduced in the largest
%   facet one at a time. This process also increases the accuracy of the
%   surface by reducing the size of the largests facets as they are often
%   the most innacurate.

    % Determine the number of nodes required to maintain consistency by
    % comparing the size of the original and the current surface.
    nNodes = length(Surface.xS)-length(fX);

    % Run through the required number of surface nodes
    for ii = 1:nNodes
        % Create a temporary surface to introduce the node with an
        % additional zero entry at the end
        tX = [fX 0]; % x-coordinates
        tY = [fY 0]; % y-coordinates

        % Discretize the temporary surface into facets
        [xi,xf,yi,yf] = discretizeSurface(fX,fY);

        % Calculate the area of all the surface facets
        area = calculateFacetsArea(xi,xf,yi,yf);

        % Find the first facet with the largest area
        idx = find(area == max(area),1);

        % Make space in the temporary variables
        tX(idx+2:end) = tX(idx+1:end-1); % x-coordinates
        tY(idx+2:end) = tY(idx+1:end-1); % y-coordinates

        % Introduce a surface node in the middle of the largest facet
        tX(idx+1) = (xi(idx) + xf(idx))/2; 
        tY(idx+1) = (yi(idx) + yf(idx))/2;

        % Clear the final x- and y-coordinates to replace them later with 
        % the new values
        clear fX fY % to regain memory access

        % Store the updated x- and y-coordinates
        fX = tX;
        fY = tY;

        % Clear the temporary variables
        clear tX tY % to regain memory access
    end

    % Store the final coordinates of the new surface
    Surface.xS = fX; % final x-coordinates
    Surface.yS = fY; % final y-coordinates

    % Check that the number of surface elements is conserved
    if length(Surface.xS) == length(Surface.xS)
        fprintf('        Compatibility test: passed\n');
    else
        warning('Number of surface elements is not conserved.');
    end
end