function kk = thermalConductivity_old(r,s)
%       Input a z cootdinate in the snow and return the conductivity at
%  that coordinate.
    
    zz = r.y; % [m] z coordinates
    
    % Make it vary linearly from 0.1-1 peak to bottom
    kk = -0.9/3/0.0254*zz + 1; % [W/m-K]

end