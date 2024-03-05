clc
clear all
close all
addpath(genpath("Boundary Conditions")); % workspace folders
addpath(genpath("Molecular Path")); % workspace folders

% Define the testing parameters
in.xS = [0 1 2];
in.zS = [0.5 1 0.5];
in.surfTemperatures = [150, 150];

SMC.material = 'ice';
SMC.FNUM = 1e16;
SMC.time = 8640;

% Assertions
tol = 1e-4;
success = 0;
failures = 0;

%% Unit Tests %%
%% Surface Data
[sData,SMC] = surfaceData(in,SMC);
% clc
% Assertions
% Surface area
trueArea = 1.1180;
if (trueArea - sData.areaFacet(1))/trueArea < tol
    success = success + 1;
else
    failures = failures + 1;
end

% Molar mass
trueMass = 2.9915e-26;
if (trueMass - sData.mass(1))/trueMass < tol
    success = success + 1;
else
    failures = failures + 1;
end

% Slope
trueSlope = 0.5;
if (trueSlope - sData.mFacet(1))/trueSlope < tol
    success = success + 1;
else
    failures = failures + 1;
end

% Y-intersect
trueY = 0.5;
if (trueY - sData.nFacet(1))/trueY < tol
    success = success + 1;
else
    failures = failures + 1;
end

% Number of particles
trueParticles = 59600;
if (trueParticles - sData.nParticles(1))/trueParticles < tol
    success = success + 1;
else
    failures = failures + 1;
end

%% Number Flux
[nFlux,Mw] = numberFlux(SMC.material,in.surfTemperatures);

% Number Flux
truenFlux = 59600;
if (truenFlux - nFlux/truenFlux) < tol
    success = success + 1;
else
    failures = failures + 1;
end

% Internal Molar Mass
truenMw = 2.9915e-26;
if (truenMw - Mw/truenMw) < tol
    success = success + 1;
else
    failures = failures + 1;
end

%% Saturation Vapor Pressure

satP = saturationPressure(T)