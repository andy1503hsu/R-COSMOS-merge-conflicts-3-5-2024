function data = monteCarloData(in,sD,SMC) %#codegen
%       Assigns a random position and launch angle to every single
%   simulated molecule bundle. The launch angles are extracted from a
%   isotropic distribution. The data is stored in a column major array.

    % Random Number Generator
    rng shuffle

    % Boltzman Constant
    kB = 1.380649e-23; % [J/K]

    % Molecular mass of the species
    Mw = molecularMass(SMC.species); % [kg/kmol]

    % Constants and preallocation
    nP = SMC.nParticles * ones(1,length(sD.xi)); % number of particles
    data.xL = zeros(1,sum(nP)); % X launch coordinates
    data.zL = zeros(1,sum(nP)); % Z launch coordinates
    data.Psi = zeros(1,sum(nP)); % projected 3D launch angles
    data.iL = zeros(1,sum(nP)); % launch facet index
    data.vL = zeros(1,sum(nP)); % launch velocities
    data.vX = zeros(1,sum(nP)); % horizontal projected velocities
    data.vZ = zeros(1,sum(nP)); % vertical projected velocities

    factor = 0; % multiplication factor used for the initial facet
    for ii = 1:length(sD.xi)

        % Index for the locations of all particles in the current facet
        index = (factor*sum(nP(1:ii-1))+1):sum(nP(1:ii));
        factor = 1; % factor used for the remaining facets

        % Constants
        dX = sD.xf(ii)-sD.xi(ii); % delta X
        dZ = sD.zf(ii)-sD.zi(ii); % delta Z

        % Initial Position of the Molecules
        data.xL(index) = sD.xi(ii) + dX.*rand(1,nP(ii)); % rng X

        if sD.mFacet(ii) == Inf % for undefined slopes
            data.zL(index) = sD.zi(ii) + dZ*rand(1,nP(ii));
        else % evaluate each X coordinate to obtain Z coordinate
            data.zL(index) = sD.mFacet(ii).*data.xL(index) + sD.nFacet(ii);
        end

        % Isotropic launch angles
        % Theta: Angle w.r.t the normal vector of the facet
        r = rand(1,nP(ii)); % randNumber for Theta
        %  cos(Theta)            sin(Theta)
        cT = (1-r).^0.5;        sT = r.^0.5;

        % Phi: Angle that rotates around the normal vector to the facet
        Phi = 2*pi*rand(1,nP(ii)); % CDF Phi
        sP = sin(Phi); % sin(Phi)

        % Psi: Projected 3D angle into the XY plane
        sTsP = sT.*sP; % constant
        N = sTsP*(dZ) + cT*(dX); % numerator
        D = sTsP*(dX) - cT*(dZ); % denominator
        data.Psi(index) = atan2d(N,D); % [deg] projected 2D launch angles

        % Launch velocities
        % Random numbers for launch velocities
        r1 = rand(1,nP(ii)); r2 = rand(1,nP(ii)); r3 = rand(1,nP(ii));

        % Maxwell distribution direct sampling
        u = (-2*log(r1)-2*log(r2).*sin(2*pi*r3).^2).^0.5; % probability
        vL = u*(kB*in.tempProfile(ii)/Mw)^0.5; % [m/s] scale
        data.vL(index) = vL; % launch velocity

        % Horizontal launch velocity in the 2D domain
        data.vX(index) = vL.*sTsP;

        % Vertical launch velocity
        data.vZ(index) = vL.*cT;

        % Launch facet index
        data.iL(index) = ii;
    end
end