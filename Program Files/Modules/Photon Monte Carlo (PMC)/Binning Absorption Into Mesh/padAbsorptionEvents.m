function absorptionEvents3d = padAbsorptionEvents(absorptionEvents, numWorkers)

    numEvents = length(absorptionEvents(:, 1));
    eventsPerWorker = ceil(numEvents / numWorkers);
    absorptionEvents3d = zeros(eventsPerWorker, 3, numWorkers);

    if mod(numEvents, numWorkers) == 0  % All workers have equal load

        for i = 1:numWorkers
            absorptionEvents3d(:, :, i) = absorptionEvents(1+(i-1)*eventsPerWorker:i*eventsPerWorker, :);
        end

    else  % Some workers will have 1 more (actual) event than other workers

        % All workers will either have minLoad or minLoad + 1 events
        minLoad = floor(numEvents / numWorkers);
        numRemainderEvents = rem(numEvents, numWorkers);
        for i = 1:numWorkers
            absorptionEvents3d(1:end-1, :, i) = absorptionEvents(1+(i-1)*minLoad:i*minLoad, :);
            if i <= numRemainderEvents % Add "remainder" events to workers
                absorptionEvents3d(end, :, i) = absorptionEvents(end-i+1, :);
            end
        end
    end

end