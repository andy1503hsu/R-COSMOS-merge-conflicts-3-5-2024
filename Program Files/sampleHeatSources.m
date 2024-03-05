function eQ = sampleHeatSources(model,data,source,dd,ss,xQ,zQ)
%       Sample heat sources at the query points (xQ,zQ) for all the day and
%   step combos specified.
    
    eQ = zeros(length(dd),length(xQ)); % allocate space
    
    % Grab the energy absorption from the specified radiation source
    if lower(source) == "light"
        eA = data.mesh.lightAbsorption;
    elseif lower(source) == "thermal"
        eA = data.mesh.thermalAbsorption;
    else
        error("Unkown radiation source specified");
    end
    
    dt = model.time.dayLength/model.time.daySteps; % [s] time step size
    for ii = 1:length(dd)
        
        % Extract nodes and elements for the new day
        ele = data.mesh.elements{dd(ii)}; % elements
        nds = data.mesh.nodes{dd(ii)}; % nodes
        
        msh = triangulation(ele,nds);
        id = pointLocation(msh,xQ',zQ');
        
        [eQ_,~,~,~] = getHeatSources(dt,ele,nds,eA{dd(ii)}(:,ss(ii)));
        
%         eQ(ii,:) = griddata(eCx,eCz,eQ_,xQ,zQ,'linear');
        eQ(ii,:) = eQ_(id);
    end
end