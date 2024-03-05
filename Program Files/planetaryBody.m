function model = planetaryBody(model)
%       Constructs a model of the planetary body where the surface is 
%   located using the properties selected by the user: planetary body mass,
%   equatorial radius, and rotation period (also defined here as the length
%   of day).

    % Mass of the planetary body
    model.body.mass = planetaryMass(model.body.mass); % [kg]

    % Equatorial radius
    model.body.radius = equatorialRadius(model.body.radius); % [m]

    % 
    model.body.rotation = rotationPeriod(model.body.rotation); % [s]

end