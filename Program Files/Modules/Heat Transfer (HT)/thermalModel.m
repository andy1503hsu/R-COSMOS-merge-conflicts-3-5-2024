function temp = thermalModel(temp,HT)
%       Create a "polygon" of the thermal domain and generate its finite
%   element mesh. Return the thermal pde model object with the mesh stored
%   inside.
%       For the "Planetary case", also returns the nodes and elements for
%   the center region of the mesh. This subset is used to bin photon
%   absorption in the PMC module.

    % Surface cordinate set
    xS = temp.xS; % x surface coordinates
    zS = temp.zS; % z surface coordinates

    % Surface coordinates
    %    Initial                      Final
    xi = xS(1:end-1);        xf = xS(2:end); % [m] x coordinates 
    zi = zS(1:end-1);        zf = zS(2:end); % [m] z coordinates

    % Area of the facets (length for 2D model)
    areaFacet = ((xf-xi).^2 + (zf-zi).^2).^0.5; % distance formula
    sN = length(areaFacet); % number of surface facets

    % Planetary large-scale
    if HT.type == "Planetary"
    
        wD = max(xS) - min(xS); % width of domain
        
        % Duplicate the granular region on each side and store the x and z
        % coordinates of the polygon vertices
        vs = [[xS-wD xS(2:end-1) xS+wD xS(end)+[wD 0] xS(1)-[0 wD wD]];...
              [zS zS(2:end-1) zS 0 0 0 0 zS(1)]];
        
        % Geometry description matrix
        gd = zeros(12,sN*3+7); % matrix to store geometry info
        gd(1,:) = 2*ones(1,sN*3+7); % 2 specifies line segments
        gd([2,4],:) = vs(:,[1:sN*3+5,sN+1,2*sN+1]); % start coordinates
        gd([3,5],:) = vs(:,[2:sN*3+6,3*sN+4,3*sN+3]); % end coordinates
        gd(6,:) = [zeros(1,sN*3+5),2,3]; % regions to the left
        gd(7,:) = [ones(1,sN),2*ones(1,sN),3*ones(1,sN),3,3,2,1,1,1,2]; % R
        
    % ARK chamber | Ocean Worlds Lab, real-world (meter-scale)
    elseif HT.type == "ArkChamber" || HT.type == "Open"
    
        % Sort the coordinates of vertices that make up the domain
        xD = [xS,xS(1)];
        zD = [zS,zS(1)];
        
        % Construct heat transfer geometry description matrix
        gd = zeros(12,sN+1); % matrix to store geometry info
        gd(1,:) = 2; % 2 specifies line segments
        gd([2,4],:) = [xD(1:end-1);zD(1:end-1)]; % start coordinates
        gd([3,5],:) = [xD(2:end)  ;zD(2:end)  ]; % end coordinates
        gd(6,:) = 0; % regions to the left
        gd(7,:) = 1; % regions to the right
    
    % Microscopic laboratory samples, ice crystals, comets
    elseif HT.type == "Free"
    
        % Geometry description matrix
        gd = zeros(12,sN); % matrix to store geometry info
        gd(1,:) = 2; % 2 specifies line segments
        gd([2,4],:) = [xi;zi]; % start coords of line segs
        gd([3,5],:) = [xf  ;zf ]; % end coords of line segs
        gd(6,:) = 1; % regions to the left
        gd(7,:) = 0; % regions to the right
    
    else
        error('Simulation type ''%s'' is not installed.',HT.type)
    end
    
    % Create thermal PDE object
    temp.tModel = createpde("thermal","transient"); % transient analysis
    
    % Assign geometry to thermal model
    geometryFromEdges(temp.tModel,gd);
    
    % Create mesh from geometry
    generateMesh(temp.tModel,... % add thermal model
                'GeometricOrder','linear',... % 3-noded linear triangles
                'Hmin',min(areaFacet)*1e-1,... % miminum edge length
                'Hmax',HT.Hmax,... % maximum edge length
                'Hgrad',HT.Hgrad); % mesh element growth rate
    
    % Assign thermal properties
    k = @(l,s) thermalConductivity(l,s,HT.rho,HT.I,HT.species);
    % Wrap heat capacity function into (location,state) definition
    Cp = @(l,s) specificHeat(s.u,HT.species);

    thermalProperties(temp.tModel,'ThermalConductivity',k,... % [W/m-K]
                             'MassDensity',HT.rho,... % [kg/m^3]
                             'SpecificHeat',Cp); % [J/kg-K]

    if HT.type == "Planetary"
        % Assign zero heat flux to the floor boundary
        thermalBC(temp.tModel,'Edge',3*sN+(2:4),'HeatFlux',0);

    elseif HT.type == "ArkChamber" || HT.type == "Open"
        % Assign a fixed temperature to the floor boundary
        thermalBC(temp.tModel,'Edge',sN+1,'Temperature',HT.walT(1)); % [K]
    end

    % Extract the nodes and elements of the mesh
    if HT.type == "Planetary" % for planetary
        
        % Extract center mesh
        [nds,ele,~] = extractCenterMesh(temp.tModel); % center region
        
        % Find the nodes and elements
        msh = triangulation(ele,nds);
        temp.meshNodes = msh.Points; % nodes
        temp.meshElements = msh.ConnectivityList; % elements

    else % for any other simulation type
        temp.meshNodes = temp.tModel.Mesh.Nodes; % nodes
        temp.meshElements = temp.tModel.Mesh.Elements; % elements
    end
end