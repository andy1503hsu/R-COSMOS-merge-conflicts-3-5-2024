function [nds,ele,nId] = extractCenterMesh(tModel)
%       Extract the nodes and elements of the planetary mesh that are
%   located in the center region of the heat transfer geometry.
    
    ndsAll = (tModel.Mesh.Nodes)'; % nodes in the entire mesh
    eleAll = (tModel.Mesh.Elements)'; % elements in the entire mesh
    
    % Find indices of nodes and elements in the center region (face 2)
    nId = findNodes(tModel.Mesh,'region','face',2); % node ids
    eId = findElements(tModel.Mesh,'region','face',2); % element ids
    
    ndOrd(nId) = 1:length(nId); % vector specifying new node order
    
    nds = ndsAll(nId,:); % the nodes in the new mesh
    ele = ndOrd(eleAll(eId,:)); % the elements in the new mesh
    
end