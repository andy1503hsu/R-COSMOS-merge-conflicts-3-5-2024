function energyInMesh = binParallelized(aEvents, nodes, elements) %#codegen

    % aEvents is either a cell array vector (thermal radiation) or a 3D
    % matrix (solar radiation)

    if iscell(aEvents) % Thermal radiation
        numIterations = length(aEvents); % Number of facets
        thermal = 1;
    else % Solar radiation
        numIterations = size(aEvents, 3);  % # of parfor workers
        thermal = 0;
    end

    energyInMeshPerIteration = zeros(size(elements, 1), numIterations);

    % Two parfor loops (within an if-else block are used because using 1
    % single parfor loop with an if-else block inside results in aEvents
    % becoming a broadcast variable, rather than a sliced variable.
    if thermal
        parfor i = 1:numIterations
            aEvents_ = aEvents{i};
            energyInMesh_ = binningQuadTree(aEvents_, nodes, elements, 1e7);
            energyInMeshPerIteration(:, i) = energyInMesh_;
        end

        % Energy in mesh for every facet of the snow surface
        energyInMesh = energyInMeshPerIteration; 
    
    else
        parfor i = 1:numIterations
            aEvents_ = squeeze(aEvents(:, :, i));
            % Last entry is 0s, relevant events includes padded 0s
            if ~any(aEvents_(end, :))
                lastEvent = find(aEvents_(:,3) > 0, 1, "last");
                if isempty(lastEvent)
                    aEvents_ = zeros(0, 3);
                else
                    lastEvent = lastEvent(1);  % Guarantees scalar for codegen
                    aEvents_ = aEvents_(1:lastEvent, :);
                end
            end
            energyInMesh_ = binningQuadTree(aEvents_, nodes, elements, 1e7);
            energyInMeshPerIteration(:, i) = energyInMesh_;
        end
        energyInMesh = sum(energyInMeshPerIteration, 2);
    end

end