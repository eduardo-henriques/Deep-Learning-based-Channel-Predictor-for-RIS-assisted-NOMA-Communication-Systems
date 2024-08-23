% Load the dataset
data = load('dataset_1.mat');
Y = data.Y;
K = size(Y, 3); % Number of time steps
SNR_values = size(Y, 4); % Number of SNR values
time_steps = 1:K; % Time step axis

% Extract and normalize the amplitude of the received signals for both users
amplitude_user1 = squeeze(Y(1, 1, :, 9, 1)).^2; % Amplitudes for user 1
amplitude_user2 = squeeze(Y(1, 2, :, 9, 1)).^2; % Amplitudes for user 2

max_amplitude = max(amplitude_user1(:)); % Maximum amplitude for normalization
normalized_amplitude_user1 = amplitude_user1 / max_amplitude;
normalized_amplitude_user2 = amplitude_user2 / max_amplitude;

% Plotting settings
savefig = 0;
% root_save = 'C:\Your\Desired\Path\';  Example path (You might want to write it yourself). Root path to save figures
linewidth = 2;
markersize = 10;
fontname = 'Times New Roman';
fontsize = 20;
legend_alg = {'User 1', 'User 2'};
location_4 = 'southeast';

colours = [0.0000 0.4470 0.7410;
           0.8500 0.3250 0.0980];

% Create the figure for plotting
figure;
set(gcf, 'Position', [0 0 800 600]);
semilogy(time_steps, normalized_amplitude_user1, 'Color', colours(1,:), 'LineWidth', linewidth);
hold on;
semilogy(time_steps, normalized_amplitude_user2, 'Color', colours(2,:), 'LineWidth', linewidth);

% Label axes
xlabel('Time (in samples)', 'FontName', fontname, 'FontSize', fontsize);
ylabel('Normalized Received Power', 'FontName', fontname, 'FontSize', fontsize);

% Customize the legend
legend(legend_alg, 'FontName', fontname, 'FontSize', fontsize, 'Location', location_4);
legend boxoff;
set(gca, 'FontName', fontname, 'FontSize', fontsize);

% Save the figure if required
if savefig == 1
    filename = fullfile(root_save, 'normalized_received_power_vs_timesteps_0db');
    saveas(gcf, filename, 'fig');
    saveas(gcf, filename, 'png');
    saveas(gcf, filename, 'epsc2');
end