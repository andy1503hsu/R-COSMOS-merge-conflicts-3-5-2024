function [facets,I] = adjustFlat(facets,flat,dispFacet,I) %#codegen
%% DESCRIPTION:
%       Function adjusts a flat max/min according to the displacement of
%   itself relative to the original surface. The flat surface should move 
%   straight up or down according to the variable "dispFacet" which 
%   contains the sign and magnitude of the displacement of the flat 
%   surface. If the sign is positive, the flat surface will move up, if it 
%   is negative, it will move down. A flat surface might also be described 
%   as a mesa surface.
%
%       Input variables:
%           facets --> array containing all the facets
%           flat --> index for the location of the flat surface
%           dispFacet --> displacement of the facet in space
%           I --> global index
%
%       Notes:
%           The meaning of the words "done" or "over" in this function is
%           that a facet is no longer in use and needs to be removed.

%           point 1 --> the INITIAL point on a facet to the LEFT of the 
%                       flat surface.
%           point 2 --> the FINAL point on a facet to the LEFT of the flat
%                       surface.

%           point 3 --> the INITIAL point on a facet to the RIGHT of the
%                       flat surface.
%           point 4 --> the FINAL point on a facet to the RIGHT of the flat
%                       surface.

%           point 5 --> the INITIAL point of the flat surface
%           point 6 --> the FINAL point of the flat surface.

                                %% Steps %%
%   1. Check if the flat surface has receded or expanded.
%       1.1 Find the point at which the surface converges and calculate the
%       initial and final points of the adjacent facets.

% -------------------------------------------------------------------------
%   2. Check if there is still a flat surface or if the flat surface is
%   gone.
%       2.1 If the flat surface is still there, adjust the adjacent facets 
%       to intersect with the flat surface and create a 'flat mountain'
%           2.1.1 Surface CONTRACTING: Check the intersection between the
%           adjacent facets and the flat surface and if it is ABOVE the
%           flat surface it means that the flat surface still exists.
%           2.1.2 Surface EXPANDING: Check the intersection between the
%           adjacent facets and the flat surface and if it is BELOW the
%           flat surface it means that the flat surface still exits.

%       2.2 If the flat surface is gone, it should be treated as a minimun
%       or a maximun. This algorith is similar for both cases.

%% Adjust the flat surface
    %%% Check if the flat surface recessed or expanded
    if dispFacet(flat-I) < 0 || dispFacet(flat-I) > 0
        %%% Find the point at which the surface converges
        % This loop moves from the flat surface to the left and to the
        % right simultaneously. It skipps the flat plate
        loc = 0; % initialize location index
        for i = 1:length(facets)/2
            loc = i; % store index in location variable
            % X coordinate of the FIRST point in the first facet
            p1X = facets(1,flat-I-loc);
            % X coordinate of the LAST point in the second facet
            p4X = facets(3,flat-I+loc);
            
            % Check if the first facet has been reached
            if p1X > p4X && p1X == facets(1,1)
                facets(:,flat-I-loc) = [];
                I = I + 1;
                return
            end
            % Check if the last facet has been reached
            if p1X > p4X && p4X == facets(3,end)
                facets(:,flat-I+loc) = [];
                I = I + 1;
                return
            end
            % If the X coordinates are in the normal position
            if p1X < p4X
                break
            end            
        end
        % Calculate the additional X and Y coordinates of the points in the
        % facets
        p1X = facets(1,flat-I-loc); p1Y = facets(2,flat-I-loc);
        p2X = facets(3,flat-I-loc); p2Y = facets(4,flat-I-loc);
        p3X = facets(1,flat-I+loc); p3Y = facets(2,flat-I+loc);
        p4X = facets(3,flat-I+loc); p4Y = facets(4,flat-I+loc);
        
        %%% Find whether the flat surface still exists, or if it is already
        %%% gone.
        % If loc equals 1, it means that the flat surface MIGHT still be 
        % there; therefore, check the intersection between adjacent facets
        if loc == 1
            x = [p1X p2X];
            y = [p1Y p2Y];
            xS = [p3X p4X];
            yS = [p3Y p4Y];
            [xSc1,ySc1] = intersectSurface2(x,y,xS,yS);
            
            %%% If there is no intersection there is still a flat surface.
            if isempty(ySc1)
                % Find the intersection of each adjacent facet with the
                % flat surface
                
                % Points 5 and 6 contain the X and Y coordinates of the
                % initial and final points of the flat surface
                p5X = facets(1,flat-I); p5Y = facets(2,flat-I);
                p6X = facets(3,flat-I); p6Y = facets(4,flat-I);
                
                % Check and adjust the first adjacent facet
                x = [p1X p2X];
                y = [p1Y p2Y];
                xS = [p5X p6X];
                yS = [p5Y p6Y];
                [xSc,ySc] = intersectSurface2(x,y,xS,yS);

                % If there is an intersection between the flat surface and
                % the facet on the LEFT
                if ~(isempty(xSc))
                    facets(3,flat-I-loc) = xSc(1,1); facets(4,flat-I-loc) = ySc(1,1); %point 2
                    facets(1,flat-I) = xSc(1,1);   facets(2,flat-I) = ySc(1,1); %point 5
                end
                
                % Check and adjust the second adjacent facet
                x = [p3X p4X];
                y = [p3Y p4Y];
                xS = [p5X p6X];
                yS = [p5Y p6Y];
                [xSc,ySc] = intersectSurface2(xS,yS,x,y);
                
                % If there is an intersection between the flat surface and
                % the facet on the RIGHT
                if ~(isempty(xSc))
                    facets(1,flat-I+loc) = xSc(1,1); facets(2,flat-I+loc) = ySc(1,1); %point 3
                    facets(3,flat-I) = xSc(1,1); facets(4,flat-I) = ySc(1,1); %point 6
                end
                
            %%% If there is an intersection between the two adjacent facets
            elseif ~(isempty(ySc1))
                %%% If the surface is contracting
                if dispFacet(flat - I) < 0
                    %%% If the intersection is above the flat surface it
                    %%% means that the flat surface is in fact still there
                    if ySc1 > facets(2,flat-I)
                        % Find the intersection of each adjacent facet with
                        % the flat surface.
                        % Calculate the points of the flat surface
                        p5X = facets(1,flat-I); p5Y = facets(2,flat-I);
                        p6X = facets(3,flat-I); p6Y = facets(4,flat-I);

                        % Check and adjust the first adjacent facet
                        x = [p1X p2X];
                        y = [p1Y p2Y];
                        xS = [p5X p6X];
                        yS = [p5Y p6Y];
                        [xSc,ySc] = intersectSurface2(x,y,xS,yS);
                        
                        % Adjust the first adjacent facet and the flat
                        % surface
                        if ~(isempty(xSc))
                            facets(3,flat-I-loc) = xSc(1,1); facets(4,flat-I-loc) = ySc(1,1); %point 2
                            facets(1,flat-I) = xSc(1,1); facets(2,flat-I) = ySc(1,1); %point 5
                        end
                        % Check and adjust the second adjacent facet
                        x = [p3X p4X];
                        y = [p3Y p4Y];
                        xS = [p5X p6X];
                        yS = [p5Y p6Y];
                        [xSc,ySc] = intersectSurface2(xS,yS,x,y);
                        
                        % Adjust the second adjacent facet and the flat
                        % surface
                        if ~(isempty(xSc))
                            facets(1,flat-I+loc) = xSc(1,1); facets(2,flat-I+loc) = ySc(1,1); %point 3
                            facets(3,flat-I) = xSc(1,1); facets(4,flat-I) = ySc(1,1); %point 6
                        end
                        %Adjust the second adjacent facets to be realistic
                        if facets(4,flat-I-loc-1) > facets(2,flat-I-loc)
                            facets(4,flat-I-loc-1) = facets(2,flat-I-loc);
                        
                        elseif facets(2,flat-I+loc+1) > facets(4,flat-I+loc)
                            facets(2,flat-I+loc+1) = facets(4,flat-I+loc);
                        end
                        
                    %%% However, if the intersection is below or at the
                    %%% same height than the flat surface, it means that 
                    %%% the flat surface is gone
                    elseif ySc1 <= facets(2,flat-I)
                        facets(3,flat-I-loc) = xSc1(1,1); facets(4,flat-I-loc) = ySc1(1,1); %point 2
                        facets(1,flat-I+loc) = xSc1(1,1); facets(2,flat-I+loc) = ySc1(1,1); %point 3
                        
                        %Adjust the second adjacent facets to be realistic
                        if facets(4,flat-I-loc-1) > facets(2,flat-I-loc)
                            facets(4,flat-I-loc-1) = facets(2,flat-I-loc);
                        
                        elseif facets(2,flat-I+loc+1) > facets(4,flat-I+loc)
                            facets(2,flat-I+loc+1) = facets(4,flat-I+loc);
                        end
                        
                        % Remove the cells that contain the facets that are
                        % no longer in use, in this case, only the flat
                        % surface
                        for k = 1:1
                            facets(:,flat-I) = [];
                            I = I + 1; % update global index
                        end
                    end
                    
                %%% If the surface is expanding
                elseif dispFacet(flat - I) > 0
                    %%% If there is an i    ntersection and it is above the
                    %%% flat surface it means that the flat surface is gone
                    if ySc1 > facets(2,flat-I)
                        facets(3,flat-I-loc) = xSc1(1,1); facets(4,flat-I-loc) = ySc1(1,1); %point 2
                        facets(1,flat-I+loc) = xSc1(1,1); facets(2,flat-I+loc) = ySc1(1,1); %point 3

                        %Adjust the second adjacent facets to be realistic
                        if facets(4,flat-I-loc-1) > facets(2,flat-I-loc)
                            facets(4,flat-I-loc-1) = facets(2,flat-I-loc);
                        
                        elseif facets(2,flat-I+loc+1) > facets(4,flat-I+loc)
                            facets(2,flat-I+loc+1) = facets(4,flat-I+loc);
                        end
                        
                        % Remove the cells that contain the facets that are
                        % no longer in use, in this case, only the flat
                        % surface
                        facets(:,flat-I) = [];
                        I = I + 1;
                        
                    %%% If there is an intersection and it is below the
                    %%% flat surface, the flat surface is still there
                    elseif ySc1 < facets(2,flat-I)
                        % Find the intersection of each adjacent facet with
                        % flat surface.
                        p5X = facets(1,flat-I); p5Y = facets(2,flat-I);
                        p6X = facets(3,flat-I); p6Y = facets(4,flat-I);
                        
                        % Check and adjust the first adjacent facet
                        x = [p1X p2X];
                        y = [p1Y p2Y];
                        xS = [p5X p6X];
                        yS = [p5Y p6Y];
                        [xSc,ySc] = intersectSurface2(x,y,xS,yS);
                        
                        % Adjust the first adjacent facet and the flat
                        % surface
                        if ~(isempty(xSc))
                            facets(3,flat-I-loc) = xSc(1,1); facets(4,flat-I-loc) = ySc(1,1); %point 2
                            facets(1,flat-I) = xSc(1,1); facets(2,flat-I) = ySc(1,1); %point 5
                        end
                        % Check and adjust the second adjacent facet
                        x = [p3X p4X];
                        y = [p3Y p4Y];
                        xS = [p5X p6X];
                        yS = [p5Y p6Y];
                        [xSc,ySc] = intersectSurface2(xS,yS,x,y);
                        
                        % Adjust the second adjacent facet and the flat
                        % surface
                        if ~(isempty(xSc))
                            facets(1,flat-I+loc) = xSc(1,1); facets(2,flat-I+loc) = ySc(1,1); %point 3
                            facets(3,flat-I) = xSc(1,1); facets(4,flat-I) = ySc(1,1); %point 6
                        end
                        % Adjust the second adjacent facets to be realistic
                        if facets(4,flat-I-loc-1) < facets(2,flat-I-loc)
                            facets(4,flat-I-loc-1) = facets(2,flat-I-loc);
                        
                        elseif facets(2,flat-I+loc+1) < facets(4,flat-I+loc)
                            facets(2,flat-I+loc+1) = facets(4,flat-I+loc);
                        end
                    end
                end
            end
         
        %%% If location is bigger than 1, it means that the flat surface is gone,
        %%% therefore, it should be treated as a maximum or minimun.
        elseif loc > 1
            % Remove the flat surface
            facets(:,flat-I) = [];
            I = I + 1; % update the global index
            
            % Remove the cells that contain the facets that are not in use
            tempIndex = 0;
            for j = 1:(loc-1)*2
                tempIndex = j;
                facets(:,flat-I-loc+2) = [];
            end
            
            % Check the new peak for an intersection
            x = [p1X p2X];
            y = [p1Y p2Y];
            xS = [p3X p4X];
            yS = [p3Y p4Y];
            [xSc,ySc] = intersectSurface2(x,y,xS,yS);

            % If there is an intersection, adjust the peak
            if ~(isempty(xSc))
                facets(3,flat-I-loc+1) = xSc(1,1);
                facets(4,flat-I-loc+1) = ySc(1,1);
                facets(1,flat-I-loc+2) = xSc(1,1);
                facets(2,flat-I-loc+2) = ySc(1,1);
            end
            I = I + tempIndex; % update the global index
        end
    end
end