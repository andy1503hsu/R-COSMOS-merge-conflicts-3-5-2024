function [eQ,eV,eCx,eCz] = getHeatSources(dt,ele,nds,eleA)
    
    nX = nds(:,1); nZ = nds(:,2); % [m] nodal x and z coord
    eX = nX(ele); % [m] x coords of each element
    eZ = nZ(ele); % [m] z coords of each element
    eCx = mean(eX,2); % [m] x coords of centroid of each element
    eCz = mean(eZ,2); % [m] z coords of centroid of each element
    
    eV = (eX(:,1).*(eZ(:,2)-eZ(:,3)) +...
          eX(:,2).*(eZ(:,3)-eZ(:,1)) +...
          eX(:,3).*(eZ(:,1)-eZ(:,2)))/2; % [m^3] volume of each element
    
    eQ = eleA./eV/dt; % [W/m^3] heat source value of each element

end