% These variables are values dependent *only* on the surface coordinates,
% and independent of the photon's location.

%       xi|xf = initial|final X coordinates of snow surface
%       zi|zf = initial|final Z coordinates of snow surface

function [x34, z34, determinant34] = cachedInfo(xi, xf, zi, zf)

    x34 = xi-xf;
    z34 = zi-zf;    
    determinant34 = xi.*zf-zi.*xf;

end