function ind = adjustCriticalFacets(fLoc,sData,Nodes,dS,ind) 
%       Adjusts the surface at the critical locations by testing if the
%   facets have displaced over the convergence point between the travel
%   path of its corresponding surface nodes.

    % Run through each critical facet
    for ii = 2:length(ind)-1
        % Check if the current facet is a critical facet.
        if fLoc(ii)
            % Retrieve the travel path of the left node (Line L1)
            %    P1(x1,y1)               P2(x2,y2)
            x1 = Nodes.xi(ii-1);    x2 = Nodes.xf(ii-1); % x-coordinates
            y1 = Nodes.yi(ii-1);    y2 = Nodes.yf(ii-1); % y-coordinates

            % Retrieve the travel path of the right node (Line L2).
            %    P3(x3,y3)               P4(x4,y4)
            x3 = Nodes.xi(ii);      x4 = Nodes.xf(ii); % x-coordinates
            y3 = Nodes.yi(ii);      y4 = Nodes.yf(ii); % y-coordinates
            
            % Find the intersection point between the travel path of both
            % surface nodes.
            [x0,y0] = intersectLines(x1,y1,x2,y2,x3,y3,x4,y4);
            
            % Retrieve the initial and final coordinates of the facet.
            x1 = sData.xi(ii);      x2 = sData.xf(ii); % x-coordinates
            y1 = sData.yi(ii);      y2 = sData.yf(ii); % y-coordinates
            
            % Calculate the minimum distance from a point in 2D space to a
            % line.
            distance = abs((x2-x1)*(y1-y0) - (x1-x0)*(y2-y1))/...
                                                sqrt((x2-x1)^2+(y2-y1)^2);
            
            % Check if the displacement of the facet exceeded the minimum
            % distance.
            if abs(dS(ii)) > distance
                ind(ii) = true; % mark the facet
            end
        end
    end
end