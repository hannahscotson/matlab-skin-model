function [model_FSR_noisy] = addNoise(numTrialsNoise, surroundingNodes)
    % Load stored model response
    load('model_FSR_normalised.mat', 'model_FSR');

    % Define letters and stiffness levels
    letters = fieldnames(surroundingNodes);
    hardness = {'Soft', 'Medium', 'Hard'};

    % Initialise a new cell array to store noisy trials
    model_FSR_noisy = cell(size(model_FSR, 1), numTrialsNoise, length(letters), length(hardness));

    % Add random noise with varying intensity
    minNoiseLevel = 0.01; % Minimum noise scaling factor
    maxNoiseLevel = 0.2;  % Maximum noise scaling factor

    % Iterate over letters (C, D, Q)
    for letterIdx = 1:length(letters)
        letter = letters{letterIdx};
        affectedNodes = surroundingNodes.(letter);  % Nodes for this letter

        % Iterate over stiffness levels (Soft, Medium, Hard)
        for stiffnessIdx = 1:length(hardness)
            for row = 1:size(model_FSR, 1)
                for trial = 1:numTrialsNoise
                    noiseLevel = minNoiseLevel + (maxNoiseLevel - minNoiseLevel) * rand(); % Random noise level

                    % Extract trial data
                    trialData = model_FSR{row, 1};

                    % Apply noise only to affected nodes
                    noisyData = trialData;
                    noisyData(affectedNodes) = noisyData(affectedNodes) + noiseLevel * randn(1, numel(affectedNodes));

                    % Store the modified data
                    model_FSR_noisy{row, trial, letterIdx, stiffnessIdx} = noisyData;
                end
            end
        end
    end

    % Save the noisy trials
    save('model_FSR_noise.mat', 'model_FSR_noisy');
end
