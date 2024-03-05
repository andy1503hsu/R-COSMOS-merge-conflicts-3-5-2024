function rho = density(r,s)
%       Input a z cootdinate in the snow and return the density at that
%  coordinate.
    
    zz = r.y; % [m] z coordinates
    
    % Make it vary linearly from 300-500 peak to bottom
    rho = -200/3/0.0254*zz + 500; % [kg/m^3]

end