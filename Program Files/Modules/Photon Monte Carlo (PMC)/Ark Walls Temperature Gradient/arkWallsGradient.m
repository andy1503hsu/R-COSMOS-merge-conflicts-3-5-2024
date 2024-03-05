%{
    Creates a temperature gradient along the ark chamber walls, from the floor to the SLI -- each wall
    is isothermal, but with enough walls, the gradient becomes fine
    enough and is functionally equivalent to having non-isothermal walls.

    Parameters:
        dDom -- The dimensions of the Ark Chamber, size is 1 x 2
            dDom(1) is the arkWidth, dDom(2) is the arkHeight
        
        Within the infoStruct:
        wallsPerSide -- Wall segments per side of the ark chamber
        wallsFromCornerToSLI -- Wall segments from top left (or top right) corner to the SLI covering
        floorTemperature -- Temperature of the floor (in Kelvin)
        SLITemperature -- Temperature of the SLI covering (in Kelvin)

        SLIWidth -- Width of SLI at center of ceiling

     Returns:
        wallTemperatures -- A vector that contains the temperatures of each wall segment
            The length of this vector is equal to 3 + 1 + 2*wallsPerSide + 2*wallsFromCornerToSLI
            (See Methodology for why this is)

     Methodology:
        See the 8.25.21 powerpoint under the "Andy's Updates" folder on Box for details!
        The floor of the Ark Chamber is treated as 3 wall segments, while
        the SLI covering is treated as 1 wall segment (this is where the 3
        + 1 comes in!)

    Last edited: 9/1/21 by Andy Hsu
%}

function [wallTemperatures, wallCoords, SLI_index] = arkWallsGradient(dDom, infoStruct, SLIWidth) %#codegen
    
    arkWidth = dDom(1);     arkHeight = dDom(2);
    dist_cornerToSLI = (arkWidth - SLIWidth) / 2; % dist from top left (or top right) corner to edge of SLI
    manhattanDist_FTSLI = arkHeight + dist_cornerToSLI; % manhattan dist from floor to (edge of) SLI

    % Make the function handle for the temperature gradient function
    tempFunction = tempGradientFunctions(infoStruct.floorTemperature, infoStruct.SLITemperature, 0, manhattanDist_FTSLI, infoStruct.tempGradientType, 1);

    % wall segment endpoints from bottom left corner to left edge of SLI (clockwise)
    wallEndpoints = [linspace(0, arkHeight, infoStruct.wallsPerSide + 1), linspace(arkHeight, arkHeight + dist_cornerToSLI, infoStruct.wallsFromCornerToSLI + 1)];
    wallEndpoints(infoStruct.wallsPerSide + 1) = []; % get rid of duplicate arkHeight value

    % temps of "left side" of ark chamber --> includes temps for left wall and left-of-SLI ceiling
    % portion
    wallTempsOneSide = averageValue(tempFunction, wallEndpoints(1:end - 1), wallEndpoints(2:end));

    wallTemperatures = [infoStruct.floorTemperature*ones(1, 2), wallTempsOneSide,...
                        infoStruct.SLITemperature, fliplr(wallTempsOneSide), infoStruct.floorTemperature];
    
    % Coordinates of vertices that make up ark walls, NOT including the two
    % coords where the ice cake meets the floor
    xCoords = [zeros(1, infoStruct.wallsPerSide),linspace(0,dist_cornerToSLI,infoStruct.wallsFromCornerToSLI + 1),...
          linspace(dist_cornerToSLI+SLIWidth,arkWidth,infoStruct.wallsFromCornerToSLI + 1),arkWidth*ones(1, infoStruct.wallsPerSide)]; % [m] x coords
    zCoords = [linspace(0, arkHeight, infoStruct.wallsPerSide + 1),arkHeight*ones(1,infoStruct.wallsFromCornerToSLI*2),...
           linspace(arkHeight, 0, infoStruct.wallsPerSide + 1)]; % [m] z coords
    
    wallCoords.xS = xCoords;
    wallCoords.zS = zCoords;
    SLI_index = infoStruct.wallsPerSide + infoStruct.wallsFromCornerToSLI + 1;
end
