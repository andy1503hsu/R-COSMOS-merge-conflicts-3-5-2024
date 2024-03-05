function msh = getArkMesh(sData, inPMC)
%       Create the mesh in which absorbed photon energy will be collected
%   from the same mesh used for the finite element heat transfer solver.
    
    xS = [sData.xi sData.xf(end)];  % surface x coords
    zS = [sData.zi sData.zf(end)];  % surface z coords
    Hmin = inPMC.Hmin; % min element length
    Hmax = inPMC.Hmax; % max element length
    Hgrad = inPMC.Hgrad;

    % Create a thermal model, from which we will extract the mesh
    tModel = createThermalModelArkChamber(xS,zS,Hmin,Hmax,Hgrad);
            
    % Extract the nodes and elements of the mesh
    nds = (tModel.Mesh.Nodes)'; % nodes
    ele = (tModel.Mesh.Elements)'; % elements
    
    msh = triangulation(ele,nds);
    
end