function sT = sampleSurfaceTemperatures(sId,soln)
%       Extract temperatures across each given surface segment and find
%   the average temperature of that surface.
    
    ndX = soln.Mesh.Nodes(1,:)'; % node x coords
    ndZ = soln.Mesh.Nodes(2,:)'; % node z coords
    ndT = soln.Temperature(:,end); % node temperatures
    
    sT = zeros(1,length(sId)); % allocate space for surface temperatures
    
    % Loop through each surface segment
    for ii = 1:length(sId)
        
        % Find the mesh nodes residing on the surface segment
        sN = findNodes(soln.Mesh,'region','Edge',sId(ii));
        
        [sX,id] = sort(ndX(sN)); % sort the nodes by increasing x coords
        sZ = ndZ(sN(id)); % corresponding z coords
        sD = ((sX(1)-sX).^2 + (sZ(1)-sZ).^2).^0.5; % dist from left node
        sT_ = ndT(sN); % corresponding temperatures
        
        % Average the temperatures for that surface segment
        sT(ii) = trapz(sD,sT_)/sD(end);
    end
    
end