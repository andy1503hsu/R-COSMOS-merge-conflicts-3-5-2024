function kk = thermalConductivityArk(r,s)
%       Input a z cootdinate in the snow and return the conductivity at
%  that coordinate.
    global Case kbulkice kfloor
    zz = r.y; % [m] z coordinates
    
    kkBulk = 0.9; % conductivity of the ice
    kkBulk = kbulkice(Case); % Notice this new line!!!!!!!!!!!!!
    kkJoint = 0.012; % the so called thermal joint from Dan's theory
    kkJoint = kfloor(Case); % Notice this new line!!!!!!!!!!!!!
    
    kk = kkBulk*ones(size(zz)); % [W/m-K] start with uniform conductivity
    
    id = zz < 0.002; % logical indices of zz values less than 1 mm
    kk(id) = kkJoint; % zz values less than 1mm are joint nodes

end