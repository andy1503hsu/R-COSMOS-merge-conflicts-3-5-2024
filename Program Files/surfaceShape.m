function surface = surfaceShape(surface)
%       Creates a sinusoidal snow field with a number of penitentes of
%   a certain width and height, wP and hP, at a given depth with a certain
%   roughness and shifted a certain number of elements to the right.

    % Preallocate surface shape distribution
    surface.x = 0;                  % [m] initial surface x coordinates
    surface.z = surface.depth(1);   % [m] initial surface z coordinates

    for ii = 1:length(surface.type) % for each surface type
        % User inputs
        t = surface.type(ii);       % shape of the surface
        h = surface.height(ii);     % [m] penitente height
        w = surface.width(ii);      % [m] penitente width
        d = surface.depth(ii);      % [m] depth below trough
        n = surface.number(ii);     % number of surface replicas
        e = surface.facets(ii);     % number of surface elements/replica
        s = surface.divisions(ii);  % number of facets/facet
        r = surface.roughness(ii);  % surface roughness
        p = surface.modifier(ii);   % modifying factor

        for jj = 1:n % for each surface replica
            % Surface x coords
            x = linspace(0,w,e+1);
    
            if t == "sinusoidal"
                % Surface z coords
                z = (h/2)*sin(x/w*2*pi-pi/2-p/e*2*pi) + ...
                        h/2 + d + [0,r*w*(rand(1,length(x)-2)-1/2),0];
        
            elseif t == "triangular"
                % Surface z coords
                z = d + h - h*abs(2*x/w-1) +...
                                  [0,r*w*(rand(1,length(x)-2)-1/2),0];
        
            elseif t == "convex"
                % Surface z coords
                z = d + h - h*abs(2*x/w-1).^(0.5+p) +...
                                  [0,r*w*(rand(1,length(x)-2)-1/2),0];
        
            elseif t == "concave"
                % Surface z coords
                z = d + h - h*abs(2*x/w-1).^(2+p) +...
                                  [0,r*w*(rand(1,length(x)-2)-1/2),0];
                
            else
                error('Surface type not installed.')
            end
    
            % Divide the facets into how many segments
            [x,z] = subdivideFacets(x,z,s);
    
            % Ensure surface fit
            x(1) = [];
            z(1) = [];
    
            % Create surface
            surface.x = [surface.x x + surface.x(end)];
            surface.z = [surface.z z];

        end
    end

    % Enforce periodic boundary conditions
    surface.z(end) = surface.z(1);
    % plot(surface.x,surface.z)
    % axis equal
end