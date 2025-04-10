%% Setup
clearvars
clc
% Generate a trial for model response data
generateModelFlg = 0;
% Generate a trial for model response data (temporal response, one combination)
generateSpecificModelFlg = 0;

% Complete classification for data response
dataClassFlg = 0;
% Complete classification for model response
modelClassFlg = 0;

% Generate animation for selected model response
plotFlg = 0;
% Generate heatmap for selected data response
dataMapFlg = 0;
% Generate heatmap for selected model response
modelMapFlg = 0;
% Generate heatmap for selected model response (temporal response, one combination)
modelTempMapFlg = 1;

%% Skin model parameters
% Indentation force
InputForce = -100;
% Input Gaussian Sigma
InputSigma = 1e-3;
% Mass value
InpM = 1e-2;

%% Arrange the masses as triangulated network across all the three layers
% Number of degree-of-freedom of the skin model.
nDoF = 3; 
% Size of the skin mesh ([rows columns])
nMeshNodes = [16 16];
% Size of skin = 0.056 x 0.056 x 0.013 cm^3
% Distance between the masses 
dMass = 0.003733333333; % 0.056/15 as there are 15 gaps
% Distance between the layers
dLayer = 0.0065; % 0.013/2 as there are three layers so two spaces

% Function generating the connectivity between the masses
[v,p,nNodes,~] = MeshTriangles(nMeshNodes,nDoF,dMass,dLayer);
nMass = numel(p(:,1));
nTriangles = numel(v(:,1));
nNodes.FirLay = nMass - (nNodes.SecLay+nNodes.TrdLay);

%% Generate model response for 9 combinations, 1 trial
% Note: Previously stored model response must have been deleted
if generateModelFlg == 1
    generateModel(nNodes, nMass, nDoF, nTriangles, nMeshNodes, InputForce, InpM, InputSigma, v, p);
end

%% Generate model response for 1 combinations, 1 trial, response over time
if generateSpecificModelFlg == 1
    generateModelForSpecificCondition(nNodes, nMass, nDoF, nTriangles, nMeshNodes, InputForce, InpM, InputSigma, v, p);
end
if modelTempMapFlg == 1
    letter = 'C';
    hardness = 'Medium';
    time_indices = [0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24];
    
    load('model_FSR_full_response.mat', 'reshaped_response', 'ts');

    % Use interpolation to get the closest indices for the time values
    time_idx = round(interp1(ts, 1:length(ts), time_indices, 'nearest'));

    modelHeatmapTemporal(reshaped_response, ts, time_idx, letter, hardness);
end
%% Define boundary nodes for each letter
boundaryNodesSet.C = [38 39 40 41 53 58 68 75 83 92 99 115 131 147 163 172 180 187 197 202 214 215 216 217];
boundaryNodesSet.D = [36 37 38 39 40 41 52 58 68 75 84 92 100 108 116 124 132 140 148 156 164 172 180 187 196 202 212 213 214 215 216 217];
boundaryNodesSet.Q = [38 39 40 41 42 46 53 59 61 68 76 83 91 93 99 109 115 125 131 141 147 157 163 173 180 188 197 203 214 215 216 217 218];

% Compute surrounding nodes for each letter
surroundingNodes = findSurroundingNodes(boundaryNodesSet, v);

% Access surrounding nodes for each letter
surrounding_C = surroundingNodes.C;
surrounding_D = surroundingNodes.D;
surrounding_Q = surroundingNodes.Q;

%% Call noise function
numTrialsNoise = 48;  % Define number of noise trials
model_FSR_noisy = addNoise(numTrialsNoise, surroundingNodes);

% Save noisy model
save('model_FSR_noisy.mat', 'model_FSR_noisy');

%% Perform Classification for model response
if modelClassFlg == 1
    numTrials = 25;
    accuracies = zeros(1, numTrials); % Preallocate accuracy array
    totalConfMat = zeros(9, 9); % Adjust size based on number of classes
    load('model_FSR_noise.mat', 'model_FSR_noisy');
    %load('model_FSR_trials.mat', 'model_FSR_trials')

    for i = 1:numTrials
        [accuracy, confMat] = modelClassification(model_FSR_noisy); % Run classification
        accuracies(i) = accuracy * 100; % Store accuracy in percentage
        totalConfMat = totalConfMat + confMat; % Sum confusion matrices across trials
        fprintf('Model Trial %d Accuracy: %.2f%%\n', i, accuracies(i));
    end

    % Compute statistics
    avgAccuracy = mean(accuracies);
    minAccuracy = min(accuracies);
    maxAccuracy = max(accuracies);

    % Compute correct and incorrect classifications per class
    correctPerClass = diag(totalConfMat); % True Positives (Correct)
    totalPerClass = sum(totalConfMat, 2); % Total instances of each class
    incorrectPerClass = totalPerClass - correctPerClass; % Incorrect classifications per class

    % Sort classes from most to least correctly classified
    [sortedCounts, sortedClassOrder] = sort(correctPerClass, 'descend');
    sortedIncorrectCounts = incorrectPerClass(sortedClassOrder); % Reorder incorrect classifications accordingly

    % Compute overall incorrect classifications
    totalCorrect = sum(correctPerClass);
    totalClassifications = sum(totalConfMat(:));
    totalIncorrect = totalClassifications - totalCorrect;

    % Print results
    fprintf('Model Classification: Average Accuracy over %d trials: %.2f%%\n', numTrials, avgAccuracy);
    fprintf('Minimum Accuracy: %.2f%%\n', minAccuracy);
    fprintf('Maximum Accuracy: %.2f%%\n', maxAccuracy);
    fprintf('Total Incorrect Classifications: %d\n', totalIncorrect);

    % Display class ranking based on correct classifications
    fprintf('Order of Correctly and Incorrectly Identified Classes (Most Correct First):\n');
    for j = 1:length(sortedClassOrder)
        fprintf('Class %d: %d correct, %d incorrect classifications\n', ...
            sortedClassOrder(j), sortedCounts(j), sortedIncorrectCounts(j));
    end
end

%% Generate heatmap for model response
if modelMapFlg == 1
    letter = 'C';
    hardness = 'Medium';
    
    load('model_FSR_noise.mat', 'model_FSR_noisy');
    first_trial_data = model_FSR_noisy{2,20}; % Select the combination and trial

    modelHeatmap(first_trial_data, letter, hardness);
end

%% Perform plotAnimation
if plotFlg == 1
    % Load the data stored for plotAnimation
    load('storedModelData.mat', 'storedData');
    
    % Select the trial and row from storedData for animation
    selectedRow = 1;
    selectedTrial = 1;
    
    % Extract the relevant data for plotAnimation
    x1 = storedData(selectedRow, selectedTrial).x1;
    ts = storedData(selectedRow, selectedTrial).ts;
    
    plotAnimation(p,nNodes,nMass,nDoF,nMeshNodes,x1,ts);
end

%% Perform Classification for experiment data
% Load data file
load('FSR_gscale_segment.mat');

if dataClassFlg == 1
    numTrials = 25; % Number of iterations
    accuracies = zeros(1, numTrials); % Preallocate accuracy array
    totalConfMat = zeros(9, 9); % Adjust size based on number of classes

    for i = 1:numTrials
        [accuracy, confMat] = dataClassification(); % Run classification
        accuracies(i) = accuracy * 100; % Store accuracy in percentage
        totalConfMat = totalConfMat + confMat; % Sum confusion matrices across trials
        fprintf('Trial %d Accuracy: %.2f%%\n', i, accuracies(i));
    end

    % Compute statistics
    avgAccuracy = mean(accuracies);
    minAccuracy = min(accuracies);
    maxAccuracy = max(accuracies);

    % Compute correct and incorrect classifications per class
    correctPerClass = diag(totalConfMat); % True Positives (Correct)
    totalPerClass = sum(totalConfMat, 2); % Total instances of each class
    incorrectPerClass = totalPerClass - correctPerClass; % Incorrect classifications per class

    % Sort classes from most to least correctly classified
    [sortedCounts, sortedClassOrder] = sort(correctPerClass, 'descend');
    sortedIncorrectCounts = incorrectPerClass(sortedClassOrder); % Reorder incorrect classifications accordingly

    % Compute overall incorrect classifications
    totalCorrect = sum(correctPerClass);
    totalClassifications = sum(totalConfMat(:));
    totalIncorrect = totalClassifications - totalCorrect;

    % Print results
    fprintf('Average Accuracy over %d trials: %.2f%%\n', numTrials, avgAccuracy);
    fprintf('Minimum Accuracy: %.2f%%\n', minAccuracy);
    fprintf('Maximum Accuracy: %.2f%%\n', maxAccuracy);
    fprintf('Total Incorrect Classifications: %d\n', totalIncorrect);

    % Display class ranking based on correct classifications
    fprintf('Order of Correctly and Incorrectly Identified Classes (Most Correct First):\n');
    for j = 1:length(sortedClassOrder)
        fprintf('Class %d: %d correct, %d incorrect classifications\n', ...
            sortedClassOrder(j), sortedCounts(j), sortedIncorrectCounts(j));
    end
end

%% Generate heatmap for experiment data
if dataMapFlg == 1
    load('FSR_gscale_segment.mat', 'FSR');

    selected_row = 33;  % Choose one row
    trial_data = FSR{selected_row, 1};  % Pick first trial for that row

    dataHeatmap(trial_data);
end
