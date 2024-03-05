function viewGeometryPhotonEngine(gD)
%       Display the geometry of described in the geometry description
%   matrix using a default MATLAB function.
    
    % Assemble geometry to view on MATLAB's built-in geometry viewer
    g = zeros(size(gD));
    g(1,:) = 2; % specify that all surfaces are lines
    g(2,:) = gD(1,:); % starting x coordinates
    g(3,:) = gD(1,:) + gD(3,:).*gD(5,:); % ending x coordinates
    g(4,:) = gD(2,:); % starting y coordinates
    g(5,:) = gD(2,:) + gD(4,:).*gD(5,:); % ending y coordinates
    g(6,:) = gD(6,:); % region to the left
    g(7,:) = gD(7,:); % region to the right
    
    pdegplot(g,'facelabels','on','edgelabels','on'); hold on;
    quiver(gD(1,:),gD(2,:),gD(3,:).*gD(5,:),gD(4,:).*gD(5,:),...
           'k','linewidth',2,'AutoScale','off','Color','b');
    
end