function dataHeatmap(trial_data)
    % Find max timestep
    [max_per_timestep, ~] = max(trial_data, [], 2);
    max_value = max(max_per_timestep);
    max_timesteps = find(max_per_timestep == max_value);
    
    % Find best timestep (highest average response)
    avg_fsr_per_timestep = arrayfun(@(t) mean(trial_data(t, :)), max_timesteps);
    [~, best_timestep_idx] = max(avg_fsr_per_timestep);
    % median_fsr_per_timestep = arrayfun(@(t) median(trial_data(t, :)), max_timesteps);
    % [~, best_timestep_idx] = max(median_fsr_per_timestep);
    % total_fsr_per_timestep = arrayfun(@(t) sum(trial_data(t, :)), max_timesteps);
    % [~, best_timestep_idx] = max(total_fsr_per_timestep);
    best_timestep = max_timesteps(best_timestep_idx);
    
    % Extract FSR data at best timestep
    heatmap_matrix = reshape(trial_data(best_timestep, :), [16, 16]);
    
    figure;
    imagesc(heatmap_matrix);
    colormap(jet);
    colorbar;
    title('Soft Q');
end


% %% Plot surfacemap
% figure;
% surf(heatmap_matrix);
% shading interp;
% colormap(jet);
% colorbar;
% title('FSR 3D Surface Plot');

