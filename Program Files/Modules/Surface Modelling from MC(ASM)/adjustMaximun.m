function [facets,I] = adjustMaximun(facets,maximun,dispFacet,I) %#codegen
%% DESCRIPTION
%       Function adjusts the surface at a maximun location. If the maximun
%   is receding then it is adjusted, if it expands, it is left intact.
%
%       Input variables:
%           facets --> array containing all the facets
%           maximun --> index for the location of the maximun
%           dispFacet --> displacement of the facet in space
%           I --> global index

                                %% Steps %%
%   1. Check if the surface is contracting at the location of the maximun
%   by taking the average of the displacement of the facets adjacent to the
%   maximun.

% -------------------------------------------------------------------------
%   2. Find the point at which the surface converges and calculate the
%   initial and final points of the adjacent facets.
%       2.1 Remove the cell that contains the facets that are no longer in
%       use to adjust the surface.

% -------------------------------------------------------------------------
%   3. Check the new peak for an intersection and adjust it accordingly.

%% Adjust the Maximun
    %%% If the average of facets adjacent to the maximun is negative,
    %%% meaning that they maximun is receding:
    if (dispFacet(maximun-I-1)+dispFacet(maximun-I))/2 < 0
        % Find which pair of facets has receded enough to create a new peak
        loc = 0; % initialize location index
        for i = 1:length(facets)/2
            loc = i; % store index in location variable
            % Initial X coordinate of the facet on the left
            p1X = facets(1,maximun-I-loc);
            % Final X coordinate of the facet on the right
            p4X = facets(3,maximun-I+loc-1);
            
            % Check if the first facet has been reached
            if p1X > p4X && p1X == facets(1,1)
                facets(:,maximun-I-loc) = [];
                I = I + 1;
                return
            end
            % Check if the last facet has been reached
            if p1X > p4X && p4X == facets(3,end)
                facets(:,maximun-I+loc-1) = [];
                I = I + 1;
                return
            end
            % If they are in the correct position, break the loop
            if p1X < p4X
                break
            end
        end
        % Calculate the X and Y coordinates of the points in the adjacent
        % facets.
        p1X = facets(1,maximun-I-loc);   p1Y = facets(2,maximun-I-loc);        
        p2X = facets(3,maximun-I-loc);   p2Y = facets(4,maximun-I-loc);
        p3X = facets(1,maximun-I+loc-1); p3Y = facets(2,maximun-I+loc-1);
        p4X = facets(3,maximun-I+loc-1); p4Y = facets(4,maximun-I+loc-1);
        
        % Remove the cells that contain the facets that are gone    
        tempIndex = 0;
        for j = 1:(loc-1)*2
            tempIndex = j;
            facets(:,maximun-I-loc+1) = [];
        end       
        
        % Check the new peak for an intersection
        x = [p1X p2X];
        y = [p1Y p2Y];
        xS = [p3X p4X];
        yS = [p3Y p4Y];
        [xSc,ySc] = intersectSurface2(x,y,xS,yS);

        % If there is an intersection, adjust the peak according to it
        if ~(isempty(xSc))
            facets(3,maximun-I-loc) = xSc(1,1);
            facets(4,maximun-I-loc) = ySc(1,1);
            facets(1,maximun-I-loc+1) = xSc(1,1);
            facets(2,maximun-I-loc+1) = ySc(1,1);
        end
        I = tempIndex + I; % update global index
    end
end