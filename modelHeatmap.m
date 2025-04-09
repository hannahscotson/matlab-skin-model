function modelHeatmap(first_trial_data, letter, hardness)
heatmap_matrix = reshape(first_trial_data, [16, 16]);
heatmap_matrix = (max(heatmap_matrix(:)) - heatmap_matrix) / (max(heatmap_matrix(:)) - min(heatmap_matrix(:))); % Normalise to [0 1]

figure;
imagesc(heatmap_matrix);
colormap(jet);
colorbar;
title(sprintf('%s %s', hardness, letter));
end


%% DISPLAY HEATMAP FOR DESIRED COMBINATION (MODEL) WITHOUT NOISE
% letter = 'C';
% hardness = 'Soft';
% load('model_FSR_data_normalised.mat', 'model_FSR'); % Load stored model response
% first_trial_data = model_FSR{8,3}; % Select the combination
% 
% heatmap_matrix = reshape(first_trial_data, [16, 16]);
% heatmap_matrix = (max(heatmap_matrix(:)) - heatmap_matrix) / (max(heatmap_matrix(:)) - min(heatmap_matrix(:))); % Normalise to [0 1]
% 
% figure;
% imagesc(heatmap_matrix);
% colormap(jet);
% colorbar;
% title(sprintf('%s %s', hardness, letter));