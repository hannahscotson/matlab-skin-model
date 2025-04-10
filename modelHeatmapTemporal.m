function modelHeatmapTemporal(reshaped_response, ts, time_idx, letter, hardness)
    % Loop through the time indices and generate heatmap for each timepoint
    for i = 1:length(time_idx)
        % Get the data for the specific timepoint using the correct index
        first_trial_data = reshaped_response(:, time_idx(i));  % Extract data for time i

        % Reshape data into a 16x16 matrix
        heatmap_matrix = reshape(first_trial_data, [16, 16]);

        % Normalise the data to [0, 1]
        heatmap_matrix = (heatmap_matrix - min(heatmap_matrix(:))) / (max(heatmap_matrix(:)) - min(heatmap_matrix(:)));

        % Create figure for heatmap
        figure;
        imagesc(heatmap_matrix);
        colormap(jet);
        colorbar;
        title(sprintf('%s %s at t = %.3fs', hardness, letter, ts(time_idx(i))));  % Add time value to title
    end
end