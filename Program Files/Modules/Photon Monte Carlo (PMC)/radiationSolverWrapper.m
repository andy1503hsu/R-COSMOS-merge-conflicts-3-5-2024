function temp = radiationSolverWrapper(temp, sData, PMC)

    fprintf('\nRunning Photon Monte Carlo: %s\n', PMC.type);
    startTime = tic();

    % Get mesh info and preallocate to-be-reused matrices
    temp = initializePMC(temp, sData, PMC);

    %% TODO: Streamline this better
    if PMC.type == "Ark Chamber"
        PMC = CONSTANTS_PMC(PMC);
    end

    maximumTries = 10;
    tryNum = 1;

    while tryNum <= maximumTries
        try
            fprintf("\n====== PMC Attempt %d of %d ======\n", ...
                    tryNum, maximumTries)

	        repartition = 0;
            temp = radiationSolver(temp, sData, PMC);
            break

        catch exception
            %{
            fprintf("\n\nERROR -- Saving variables for reproducibility...\n")
            folderName = "Variables during PMC Errors";
            [~, ~] = mkdir(folderName);
            filename = folderName + "/" + string(datetime()) + ".mat";
            save(filename)
            %}
            switch exception.identifier
                case 'MATLAB:badsubscript'

                % Array indices must be positive integers or logical values
                % Cause is unknown and infeasible to reproduce
                % Rerun simulation with different seed
                    fprintf("\n\nARRAY INDICES ERROR -- ")
                    fprintf("Reattempting PMC with different seed...\n\n")
                    seed = randi(2^20);
                    rng(seed)

                case 'EMLRT:runTime:NumelOverflow'

                % Max array size for MATLAB compiled-code is exceeded
                % https://www.mathworks.com/help/coder/ug/array-size-restrictions-for-code-generation.html
                % Occurs with small grain radii during solar radiation
                % Retry by partioning solar photons into smaller chunks
                    fprintf("\n\nMAX ARRAY SIZE ERROR -- ")
                    repartition = 1;

                otherwise
                    % Error with parfor due to machine running out of RAM
                    % Various error messages due to this, so otherwise 
                    % catch-all is used
                    fprintf("\n\nMACHINE RAN OUT OF RAM ERROR -- ")
		            repartition = 1;
            end
            
            if repartition
                fprintf("Reattempting PMC with more partitioning...\n\n")
                    
                % All potential partition values
                partitions = [1 2 5 10 25 50 100 250 500];
                currentIndex = find(partitions == temp.numPartitions, 1);
                if currentIndex < length(partitions)
                    newPartition = partitions(currentIndex + 1);
                    temp.numPartitions = newPartition;

                    % Recompile solar and thermal radiation solver
                    % To account for new sizing
                    refreshCompilation(temp, PMC)
                end % Repeat max partition if it comes to that
            end

            tryNum = tryNum + 1;
        end  % End of try-catch block
    end  % End of while loop
    
    if tryNum == maximumTries + 1
        fprintf("\nPMC failed after %d consecutive tries\n", maximumTries)
        fprintf("Simulation terminating...\n")
        error("PMC failed after %d consecutive tries\n", maximumTries)
    end
    fprintf("%30s %f s\n\n","Completion time:", toc(startTime));
end
