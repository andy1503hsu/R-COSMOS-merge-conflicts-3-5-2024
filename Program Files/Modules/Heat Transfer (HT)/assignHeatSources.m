function tModel = assignHeatSources(dt,ele,nds,eAL,eAB,tModel,HT)
%       Compute the internal heat sources from the absorption distributions
%   and assign them to the snow domain.

    if isempty(eAL) && isempty(eAB) % if no radiation was absorbed
        return; % no heat sources to assign
        
    elseif isempty(eAL) % if no light energy radiation absorbed
        eAL = zeros(size(ele(:,1)));
        
    elseif isempty(eAB) % if no black body radiation was absorbed
        eAB = zeros(size(ele(:,1)));
        
    end

    [eQ,~,eCx,eCz] = getHeatSources(dt,ele,nds,eAL+eAB);
    
    if HT.type == "Planetary"
        % Duplicate the heat sources to mimic periodicity
        wD = max(nds(:,1)) - min(nds(:,1)); % mesh width
        eCx = [eCx-wD;eCx;eCx+wD];
        eCz = [eCz;eCz;eCz];
        eQ = [eQ;eQ;eQ];
        
        % Get rid of overlapping nodes
        [~,idx,~] = unique([eCx,eCz],'rows');
        eCx = eCx(idx);
        eCz = eCz(idx);
        eQ = eQ(idx);
    end

    fQ = scatteredInterpolant(eCx,eCz,eQ); % heat source interp function
    
    % Assign the heat source interp function
    internalHeatSource(tModel,@(r,~) fQ(r.x,r.y));
    
%     scatter(eCx(:),eCz(:),40,eQ(:),'filled')

end