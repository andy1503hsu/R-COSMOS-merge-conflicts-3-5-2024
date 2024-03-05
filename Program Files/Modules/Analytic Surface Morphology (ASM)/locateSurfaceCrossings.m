function ind = locateSurfaceCrossings(midX,midY,ind) %#codegen
%       Finds any intersection point between non-adjacent surface facets
%   and marks the facets in between for elimination.
    
    % Discretize the surface nodes into facet coordinates
    [xi,xf,yi,yf] = discretizeSurface(midX,midY);
    
    % Run through each facet but the last 2 since those are going to be
    % checked in the nested loop
    for ii = 1:length(xi)-2
        % Retrieve the initial and final coordinates of the primary facet.
        %    Initial             Final
        x1 = xi(ii);        x2 = xf(ii); % x-coordinates
        y1 = yi(ii);        y2 = yf(ii); % y-coordinates
        
        % For each facet, run from that 
        for jj = ii+2:length(xi)
            % Retrieve the initial and final coordinates of the secondary
            % facet.
            %    Initial                 Final
            x3 = xi(jj);            x4 = xf(jj);
            y3 = yi(jj);            y4 = yf(jj);
            
            % Find the intersection between the two lines
            [pX,pY] = intersectLines(x1,y1,x2,y2,x3,y3,x4,y4);
            
            % Check if the intersection coordinates lie on both the
            % selected facets:
            %   Condition for the primary facet
            if (x2-pX)*(x1-pX) <= 0 && (y2-pY)*(y1-pY) <= 0
                %   Condition for the secondary facet
                if (x3-pX)*(x4-pX) <= 0 && (y3-pY)*(y4-pY) <= 0
                    ind(ii+2:jj) = true; % mark all the in-between facets
                end
            end
        end
    end
end