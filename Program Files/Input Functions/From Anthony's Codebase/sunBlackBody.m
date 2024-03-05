%% From's Anthony's Codebase, will be used to verify figures in Carreon paper
% For original file, go to:
% Box\Research\Code\Anthony's Codebase\For Andy\Carreon Program
% Files\Europa Model Inputs\sunIrradianceEuropa

function irr = sunBlackBody(wL,wU,dist)
%       Return the spectral irradiance of a black body at the given
%   temperature in the given wavelength bounds at the given distance dist (in
%   astronomical units). Sample nP points for the irradiance.

    % Physical constants
    c = 2.99792458e8; % [m/s] speed of light
    h = 6.62607015e-34; % [J-s] Plank constant
    k = 1.380649e-23; % [J/K] Boltzmann constant
    s = 2*pi^5*k^4/15/c^2/h^3; % [W/m^2/K^4] Stefan-Boltzmann constant

    % Relevant variables
    T = 5778; % [K] sun's surface temperature
    rS = 6.95e8; % [m] sun's radius
    rE = 1.496e11; % [m] distance from sun to Earth (equal to 1 AU)
    rB = rE*dist; % [m] Convert inputted distance to meters
    nP = 1e5; % number of data ponts
    
    ww = linspace(wL,wU,nP);

    % Plank's Law of Irradiance
    ir = 2*pi*h*(c^2)./(ww.^5)./(exp(h*c./ww/k./T)-1)*(rS/rB)^2;
    
    irr = [ww;ir]';
    
    % Compute the approximate radiant emittance and find error
    qAprox = trapz(irr(:,1),irr(:,2));
    qExact = s*T^4*(rS/rB)^2;
    
    fprintf("Black body approximation: %f %%\n",qAprox/qExact*100);
end