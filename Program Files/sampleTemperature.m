function TQ = sampleTemperature(data,dQ,sQ,xQ,zQ)
%       Sample temperatures at the query points (xQ,zQ) for all the days
%   and steps specified.

    TQ = zeros(length(dQ),length(xQ));
    
    for ii = 1:length(dQ)
        nx = data.mesh.nodes{dQ(ii)}(:,1); % [m] node z coords
        nz = data.mesh.nodes{dQ(ii)}(:,2); % [m] node x coords
        nT = data.mesh.temperatures{dQ(ii)}(:,sQ(ii)); % [K] temps
        
        TF = scatteredInterpolant(nx,nz,nT,'linear');
        TQ(ii,:) = TF(xQ,zQ);
    end
    
    
end