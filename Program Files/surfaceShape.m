function surface = surfaceShape(surface)
%       Creates a sinusoidal snow field with a number of penitentes of
%   a certain width and height, wP and hP, at a given depth with a certain
%   roughness and shifted a certain number of elements to the right.

    % User inputs
    h = surface.height;     % [m] penitente height
    w = surface.width;      % [m] penitente width
    d = surface.depth;      % [m] depth below trough
    n = surface.number;     % number of penitentes
    e = surface.facets;     % number of surface elements per penitente
    r = surface.roughness;  % surface roughness
    p = surface.modifier;      % shift the field these many elements to the right

    % Surface x coords
    surface.x = linspace(0,w*n,n*e+1);

    if surface.type == "sinusoidal"
        % Surface z coords
        surface.z = (h/2)*sin(surface.x/w*2*pi-pi/2-p/e*2*pi) + ...
                h/2 + d + [0,r*w*(rand(1,length(surface.x)-2)-1/2),0];

    elseif surface.type == "triangular"
        % Surface z coords
        surface.z = d + h - h*abs(2*surface.x/w-1) +...
                          [0,r*w*(rand(1,length(surface.x)-2)-1/2),0];

    elseif surface.type == "convex"

        % Surface z coords
        surface.z = d + h - h*abs(2*surface.x/w-1).^(0.5+p) +...
                          [0,r*w*(rand(1,length(surface.x)-2)-1/2),0];

    elseif surface.type == "concave"
        % Surface z coords
        surface.z = d + h - h*abs(2*surface.x/w-1).^(2+p) +...
                          [0,r*w*(rand(1,length(surface.x)-2)-1/2),0];
    else
        error('Surface type not installed.')
    end
end