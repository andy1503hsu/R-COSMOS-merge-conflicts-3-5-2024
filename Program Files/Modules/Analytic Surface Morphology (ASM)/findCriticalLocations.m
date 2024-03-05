function [fLoc,nLoc] = findCriticalLocations(sData,dS,Opt) %#codegen
%       Finds the critical locations along the original surface in order to
%   determine where to adjust the surface. The index of the locations
%   corresponds to the index of the facet. Maxima and minima are
%   indistinguishable from each other since their adjustment is going to be
%   determined by the sign of the surface displacement. This function also
%   returns the index of the nodes that are neither a maximum or a minimum.
    
    % Preallocate the critical locations and rotation angle arrays
    nLoc = false(1,length(sData.xi)-1); % critical locations
    fLoc = false(1,length(sData.xi)); % critical facets
    
    % Find the locations of the critical facets
    for ii = 1:length(sData.xi)-1
        % Extract the vector in the direction of the primary facet
        xP = sData.xf(ii)-sData.xi(ii);
        yP = sData.yf(ii)-sData.yi(ii);
        
        % Calculate the angle between the primary facet and the x-axis
        theta = -atan2d(yP,xP);
        
        % Extract the vector in the direction of the adjacent facet
        xA = sData.xf(ii+1)-sData.xi(ii+1);
        yA = sData.yf(ii+1)-sData.yi(ii+1);

        % Determine the angle between the adjacent facet and the x-axis
        v = [xA;
             yA];
        
        % Create a rotation matrix for a theta rotation around the z-axis
        Rz = [cosd(theta) -sind(theta);
              sind(theta)  cosd(theta)];
        
        % Rotate the adjacent facet by theta
        Rv = Rz*v;
        
        % Calculate the angle between the rotated adjacent facet by theta
        % and the x-axis
        phi = round(atan2d(Rv(2),Rv(1)),Opt.res);
        
        % Determine the maxima and minima
        % If the node is a 'maximum', only receeding facets are checked
        if phi < 0 % maximum condition.
            % Check the displacement of the left adjacent facet
            if dS(ii) < 0 % if the surface is receding
                fLoc(ii) = true; % mark the facet
            end
            
            % Check the displacement of the right adjacent facet
            if dS(ii+1) < 0 % if the surface is receding
                fLoc(ii+1) = true; % mark the facet
            end

        % If the node is a 'minimum', only expanding facets are checked
        elseif phi > 0 % minimum condition
            % Check the displacement of the left adjacent facet
            if dS(ii) > 0 % if the surface is expanding
                fLoc(ii) = true; % mark the facet
            end
            
            % Check the displacement of the right adjacent facet
            if dS(ii+1) > 0 % if the surface is expanding
                fLoc(ii+1) = true; % mark the facet
            end
            
        % Else the node is neither a maximum or a minimum
        else
            nLoc(ii) = true; % mark the node
        end
        
        % Check if none of the corresponding facets are maximum or minimum
        if ~(fLoc(ii) && fLoc(ii+1)) % if not
            nLoc(ii) = true; % mark the node
        end
    end % check the next surface node
end