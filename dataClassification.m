function [accuracy, confMat] = dataClassification()
    % Load data file
    dataFile = 'FSR_gscale_segment.mat';
    if exist(dataFile, 'file') == 2
        load(dataFile, 'FSR');
    else
        error('Data file not found: %s', dataFile);
    end

    % Define selected row indices
    selected_rows = [7, 8, 9, 19, 20, 21, 31, 32, 33];
    num_rows = length(selected_rows);
    num_trials = 48; % Each row has 48 trials

    % Initialize storage
    data = [];
    labels = [];
    label_mapping = 1:num_rows;

    % Loops over each row, extracts trial data, converts it from cell to matrix
    for i = 1:num_rows
        row_idx = selected_rows(i);
        for trial = 1:num_trials
            trial_data = cell2mat(FSR(row_idx, trial));

            % Find the timestep where the maximum first occurs
            [max_per_timestep, ~] = max(trial_data, [], 2);  % Get max per timestep
            max_value = max(max_per_timestep);  % Find overall max value

            % Find all timesteps where the maximum value occurs (may be multiple)
            max_timesteps = find(max_per_timestep == max_value); 

            % Initialize a variable to track the timestep with the highest average FSR
            avg_fsr_per_timestep = zeros(length(max_timesteps), 1);
            
            % Calculate the average FSR value across all 256 features for each of those timesteps
            for j = 1:length(max_timesteps)
                t = max_timesteps(j);
                avg_fsr_per_timestep(j) = mean(trial_data(t, :));  % Average across all 256 features
            end

            % Find the timestep with the highest average FSR value
            [~, best_timestep_idx] = max(avg_fsr_per_timestep);
            best_timestep = max_timesteps(best_timestep_idx);  % Best timestep with highest avg FSR value

            % Extract the FSR data (feature vector) at the best timestep
            feature_vector = squeeze(trial_data(best_timestep, :)); % Ensure 1x256

            % Check size of [no. of timesteps x 256]
            [t, num_features] = size(trial_data);
            % Ensure it has 256 features per timestep
            if num_features ~= 256
                warning('Skipping row %d, trial %d -> Unexpected feature size: %d (expected 256)', row_idx, trial, num_features);
                continue;
            end

            % Print size of trial matrix, selected timestep and max value
            %fprintf('Row %d, Trial %d -> Size: [%d, %d], Selected Timestep: %d, Max Value: %.4f\n', ...
            %    row_idx, trial, t, num_features, best_timestep, max_value);

            % Store the extracted features and labels
            data = [data; feature_vector]; % Append as new row
            labels = [labels; label_mapping(i)]; % Append label
        end
    end

    % Manually split data into train and test
    training = 0.8;
    numSamples = size(data, 1);
    numTrain = round(training * numSamples);

    % Randomly shuffle the data indices
    indices = randperm(numSamples);

    % Split into training and test sets
    trainIdx = indices(1:numTrain);
    testIdx = indices(numTrain+1:end);

    % Train and Test Data
    X_train = data(trainIdx, :);
    Y_train = labels(trainIdx);
    X_test = data(testIdx, :);
    Y_test = labels(testIdx);

    % Standardize the data manually
    X_train = normalize(X_train); % Standardize training data
    X_test = normalize(X_test);   % Standardize test data

    % Train SVM classifier using fitcecoc for multi-class classification
    SVMModel = fitcecoc(X_train, Y_train, 'Learners', 'Linear', 'Coding', 'onevsall');

    % Predict on test set
    Y_pred = predict(SVMModel, X_test);

    % Training/testing split
    fprintf('Training to Testing Split: %.2f : %.2f\n', training, (1-training));

    % Evaluate accuracy
    accuracy = sum(Y_pred == Y_test) / length(Y_test);
    fprintf('Data Classification Accuracy: %.2f%%\n', accuracy * 100);

    % Confusion matrix
    confMat = confusionmat(Y_test, Y_pred);
    disp('Data Confusion Matrix:');
    disp(confMat);
end
