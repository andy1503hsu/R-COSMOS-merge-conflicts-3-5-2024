function [xS,yS] = boundaryConditions(Surface,Facets,nX,nY,Opt) %#codegen
%       Applies the boundary conditions according to the simulation type.
%   For periodic simulations, all the facets with coordinates outside the
%   periodic boundaries are eliminated and the end facets are cropped to
%   the boundary.

    % For periodic boundary conditions
    if Opt.type == "Planetary"
        
        % Create a line for the left boundary
        %    Initial                     Final
        x1 = Surface.xS(1);         x2 = Surface.xS(1); % x-coordinate
        y1 = rand;                  y2 = rand; % y-coordinate
        
        % Discretize the surface nodes into facets
        [xi,xf,yi,yf] = discretizeSurface(nX,nY);
        
        % Create an index to mark facets for elimination
        idx = false(1,length(xi));
        
        % Mark all the facets outside the periodic boundaries
        idx(xi <= Surface.xS(1) & xf <= Surface.xS(1)) = true; % left
        idx(xi >= Surface.xS(end) & xf >= Surface.xS(end)) = true; % right
        
        % Eliminate the marked facets
        % Initial             Final
        xi(idx) = [];       yi(idx) = []; % x-coordinates
        xf(idx) = [];       yf(idx) = []; % y-coordinates
        
        % Retrieve the initial and final coordinates of the first facet
        %    Initial                 Final
        x3 = Facets.xi(1);      x4 = Facets.xf(1); % x-coordinate
        y3 = Facets.yi(1);      y4 = Facets.yf(1); % y-coordinate

        % Determine the intersection point between the boundary and the
        % first facet.
        [pX,pY] = intersectLines(x1,y1,x2,y2,x3,y3,x4,y4);
        

        % Determine if the intersection point lies on the first facet
        if (x3-pX)*(x4-pX) <= 0 && (y3-pY)*(y4-pY) <= 0
            % Adjust the coordinates of the facet to intersect the boundary
            Facets.xi(1) = x1; % x-coordinate
            Facets.yi(1) = pY; % y-coordinate
        end

        % Create a line for the left boundary
        %    Initial                     Final
        x1 = Surface.xS(end);       x2 = Surface.xS(end); % x-coordinate
        y1 = rand;                  y2 = rand; % y-coordinate

        % Retrieve the initial and final coordinates of the last facet
        %    Initial                     Final
        x3 = Facets.xi(end);        x4 = Facets.xf(end); % x-coordinate
        y3 = Facets.yi(end);        y4 = Facets.yf(end); % y-coordinate

        % Determine the intersection point between the boundary and the
        % last facet.
        [pX,pY] = intersectLines(x1,y1,x2,y2,x3,y3,x4,y4);

        % Determine if the intersection point lies on the last facet
        if (x3-pX)*(x4-pX) <= 0 && (y3-pY)*(y4-pY) <= 0
            Facets.xf(end) = x1; % x-coordinate
            Facets.yf(end) = pY; % y-coordinate
        end

        % Average the height of the last two facets
        avgY = (Facets.yi(1) + Facets.yf(end))/2;
        
        % Construct the complete final array including boundary conditions
        xS = [Surface.xS(1) xi xf(end) Surface.xS(end)]; % x-coordinates
        yS = [avgY yi yf(end) avgY]; % y-coordinates
    end
end