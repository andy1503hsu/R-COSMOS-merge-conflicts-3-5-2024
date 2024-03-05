function surface = surfaceProtrusion
%       Creates a sigular, sinusoidal-shaped depression at the center of a
%   flat snow field.
    
    % Parameters
    dD = -1*0.0254; % depth of depression
    wD = 0.0254; % width of depression
    nD = 20; % surface segments in the depression
    dF = 2*0.0254; % depth of flat portions
    wF = 0.5*0.0254; % width of flat portions (per side)
    nF = 2; % number of surface segments in flat portions (per side)
    
    % Construct each portion of the depression field
    xL = linspace(0,wF,nF+1); % x coords of left flat portion
    xD = linspace(wF,wF+wD,nD+1); % x coords of depression
    xR = linspace(wF+wD,2*wF+wD,nF+1); % x coords of right flat portion
    
    zL = dF*ones(size(xL)); % z coords of left flat portion
    zD = (cos(2*pi*linspace(0,1,nD+1))-1)*dD/2+dF; % z coords of depression
    zR = dF*ones(size(xR)); % z coords of right flat portion
    
    % Put it all together
    surface.x = [xL(1:end-1),xD,xR(2:end)];
    surface.z = [zL(1:end-1),zD,zR(2:end)];
    
    
%     plot(surface.x,surface.z); axis equal;
end