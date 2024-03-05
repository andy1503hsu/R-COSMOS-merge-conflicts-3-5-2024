function [zenEdges, aziEdges, zenIds, aziIds] = binZenAzi(eData, numZenBins, numAziBins)

    zenEdges = linspace(0, pi/2, numZenBins+1);
    aziEdges = linspace(0, 2*pi, numAziBins+1);
    
    [~, ~, ~, zenIds, aziIds] = histcounts2(eData(:, 3), eData(:, 4),...
                             zenEdges, aziEdges);
end