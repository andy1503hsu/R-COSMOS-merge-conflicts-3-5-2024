function rr = wallReflectivityArk(ww)
%       Return the reflectivity of the Ark chamber walls from the given
%   wavelength. Reflectivity values are based off Marchall et. al.
%   "Characterization of the Reflectivity of Various Black Materials"
%
%   INPUT: ww - [m] wavelength
%   OUTPUT: rr - [-] reflectivity

    filename = "reflectivity_anodized_6061_aluminum.mat";
    load(filename,'data');
    % yields a reflectivity of 0.2ish %
    rr = interp1(data(:,1),data(:,2),ww,'linear',data(end,2));
    
end