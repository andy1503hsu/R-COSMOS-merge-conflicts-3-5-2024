function Facets = resizeFacets(Facets,Nodes,nLoc,dS) %#codegen
%       Extends the facets to reach the path of the nodes that do not 
%   contain a minimum or a maximum.
    
    % Run through all the facets but the first and last since those will be
    % adjusted through boundary conditions 
    for ii = 2:length(Facets.xi)-1
        
        % If the facet did not move
        if dS(ii) == 0
            continue % to the next facet
        end
        
        % Retrieve the initial and final coordinates of the facet. Here
        % the facet represents line L1
        %    P1(x1,y1)                   P2(x2,y2)
        x1 = Facets.xi(ii);         x2 = Facets.xf(ii); % x-coordinate
        y1 = Facets.yi(ii);         y2 = Facets.yf(ii); % y-coordinate
        
        
        % Line L2 is represented by the travel path of the surface nodes
        % Check if the left node is not a maximum or minimum
        if nLoc(ii-1)
            
            % Retrieve the travel path of the surface node which represents
            % line L2
            %    P3(x3,y3)              P4(x4,y4)
            x3 = Nodes.xi(ii-1);   x4 = Nodes.xf(ii-1); % x-coordinate
            y3 = Nodes.yi(ii-1);   y4 = Nodes.yf(ii-1); % y-coordinate

            % Find the intersection between the two lines
            [pX,pY] = intersectLines(x1,y1,x2,y2,x3,y3,x4,y4);

            % Replace the initial coordinates of the facet by the
            % intersection point.
            Facets.xi(ii) = pX; % x-coordinate
            Facets.yi(ii) = pY; % y-coordinate
        end
        
        % Check if the right node is not a maximum or minimum
        if nLoc(ii)
            % Retrieve the travel path of the surface node which represents
            % line L2
            %    P3(x3,y3)              P4(x4,y4)
            x3 = Nodes.xi(ii);     x4 = Nodes.xf(ii); % x-coordinate
            y3 = Nodes.yi(ii);     y4 = Nodes.yf(ii); % y-coordinate

            % Find the intersection between the two lines
            [pX,pY] = intersectLines(x1,y1,x2,y2,x3,y3,x4,y4);
            
            % Replace the final coordinates of the facet by the
            % intersection point.
            Facets.xf(ii) = pX; % x-coordinate
            Facets.yf(ii) = pY; % y-coordinate
        end
    end
%     for kk = 1:length(Facets.xi)
%         hold on
%         plot([Facets.xi(kk) Facets.xf(kk)],[Facets.yi(kk) Facets.yf(kk)],'g-');
%     end
end