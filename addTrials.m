function [model_FSR_trials] = addTrials(numTrialsNoise)
    % Load stored model response
    load('model_FSR_normalised.mat', 'model_FSR');

    % Initialise a new cell array to store new trials
    model_FSR_trials = cell(size(model_FSR, 1), numTrialsNoise); 

    % Add random noise with varying intensity
    minNoiseLevel = 0; % Minimum noise scaling factor
    maxNoiseLevel = 0;  % Maximum noise scaling factor

    for row = 1:size(model_FSR, 1)
        for trial = 1:numTrialsNoise
            noiseLevel = minNoiseLevel + (maxNoiseLevel - minNoiseLevel) * rand(); % Random noise level
            
            % Extract the third column data
            thirdColumnData = model_FSR{row, 1};  
            
            % Add noise to the third column data
            model_FSR_trials{row, trial} = thirdColumnData + noiseLevel * randn(size(thirdColumnData));
        end
    end

    % Save the noisy trials
    save('model_FSR_trials.mat', 'model_FSR_trials');
end
