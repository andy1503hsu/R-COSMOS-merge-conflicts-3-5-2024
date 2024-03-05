function s = surfaceSinusoidalPenitentes
%       Creates a sinusoidal penitente field with a number of penitentes of
%   a certain width and height, wP and hP, at a given depth with a certain
%   roughness and shifted a certain number of elements to the right.

    % User inputs
    h = 0.0254;	% penitente height
    w = 0.0254;	% penitente width
    d = 2*0.0254;	% depth below trough
    n = 1;	% number of penitentes
    e = 40;	% number of surface elements per penitente
    r = 0;	% surface roughness
    p = 0;	% shift the field these many elements to the right

    % Surface x coords
    s.x = linspace(0,w*n,n*e+1);
    
    % Surface z coords
    s.z = (h/2*[1,1-r*rand(1,length(s.x)-2),1]).*...
                sin(s.x/w*2*pi-pi/2-p/e*2*pi) + h/2 + d;
            
%     plot(s.x,s.z); axis equal;
end