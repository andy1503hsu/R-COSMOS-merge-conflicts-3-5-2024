function tModel = assignTemperatures(tModel,iT,HT)
%       Assign initial temperatures in the domain and duplicate the
%   temperatures on the left and right sides for periodicity.

    if size(iT,2) ~= 3 % if the simulation is only beginning
        
        thermalIC(tModel,iT); % assign initial temperature
        return;
    end

    nX = iT(:,1); % [m] x coords of nodes
    nZ = iT(:,2); % [m] z coords of nodes
    nT = iT(:,3); % [K] temperatures

    if HT.type == "Planetary"
        % Duplicate the temperatures for periodicity
        wD = max(nX) - min(nX); % [m] mesh width
        nX = [nX-wD;nX;nX+wD];
        nZ = [nZ;nZ;nZ];
        nT = [iT(:,3);iT(:,3);iT(:,3)];
    
        % Get rid of duplicates
        [~,idx,~] = unique([nX,nZ],'rows');
        nX = nX(idx);
        nZ = nZ(idx);
        nT = nT(idx);
    end
    
    fT = scatteredInterpolant(nX,nZ,nT); % interpolant function
    thermalIC(tModel,@(r) fT(r.x,r.y)); % assign interpolant function
end