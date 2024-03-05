function [xS,yS,dS,ind] = prepareSurface(Surface,dispFacet,opt) %# codegen
%       Prepares the surface according to the simulation type selected.
%   For periodic boundary conditions, the surface is copied and added at
%   each side of the original surface. For free boundary conditions, the
%   surface is sandwiched between two copies of the original surface. For
%   open boundary conditions, two planar segments are added at the end
%   points of the surface that represent the floor.

    % For periodic boundary conditions
    if opt.type == "Planetary"
        % Preallocate new surface coordinates and displacement array
        xS = zeros(1,3*length(Surface.xS)-2); % x-coordinates
        yS = zeros(1,3*length(Surface.xS)-2); % y-coordinates
        dS = zeros(1,3*length(dispFacet)); % facets displacement
        
        % Calculate correction factor
        factor = Surface.xS(end) - Surface.xS(1);
        
        % Create the x-coordinates array
        xS(1,:) = [Surface.xS(1:end-1) - factor... % left side copy
                   Surface.xS... % original surface x-coordinates
                   Surface.xS(2:end) + factor]; % right side copy
        
        % Add the corresponding y-coordinates
        yS(1,:) = [Surface.yS(1:end-1)... % left side copy
                   Surface.yS... % original surface y-coordinates
                   Surface.yS(2:end)]; % right side copy

        % Add the corresponding surface displacement
        dS(1,:) = [dispFacet... % left side copy
                   dispFacet... % original surface displacement
                   dispFacet]; % right side copy
          
        % Preallocate a logical index to identify facets that require 
        % elimination
        ind = false(1,length(dS)); % facet elimination index
        
        % Mark the newly created facets
        ind(1:length(Surface.xS)-1) = true; % left side copy
        ind(2*length(Surface.xS)-1:end) = true; % right side copy
        
    elseif opt.type == "free"
        
        
        
    elseif opt.type == "open"
        % Preallocate new surface coordinates and displacement array
        xS = zeros(1,length(Surface.xS)+2); % x-coordinates
        yS = zeros(1,length(Surface.xS)+2); % y-coordinates
        dS = zeros(1,length(dispFacet)+2); % facets displacement
        
        % Create the x-coordinates array
        xS(1,:) = [Surface.xS(end) Surface.xS Surface.xS(1)]; % right side copy
               
        % Add the corresponding y-coordinates
        yS(1,:) = [Surface.yS(end) Surface.yS Surface.yS(1)]; % right side copy

        % Add the corresponding surface displacement
        dS(1,:) = [0 dispFacet 0]; % right side copy
        
        % Preallocate a logical index to identify facets that require 
        % elimination
        ind = false(1,length(dS)); % facet elimination index
        
        % Mark the newly created facets
        ind(1) = true; % left side copy
        ind(end) = true; % right side copy
    end
end