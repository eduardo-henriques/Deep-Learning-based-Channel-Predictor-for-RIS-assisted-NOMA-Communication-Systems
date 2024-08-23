% Load the data from the .mat files
ber_pooling_after = load('ber_epochs_100.mat');
ber_pooling_before = load('ber_old.mat');

% Extract the BER values
ber_after = ber_pooling_after.BER; % Adjust the field if necessary
ber_before = ber_pooling_before.BER; % Adjust the field if necessary

% Define the SNR values
SNR = -10:5:40; % Adjust if necessary

% Define the rows for User 1 and User 2
rows_to_plot = [3, 7];
ber_after_selected = ber_after(rows_to_plot, :);
ber_before_selected = ber_before(rows_to_plot, :);

% Define colors for each user
colors = [0.0000 0.4470 0.7410;  % Blue for User 1
          0.8500 0.3250 0.0980]; % Red for User 2

% Define markers for after and before scenarios
markers_after = {'-o', '-o'}; % Circle for after
markers_before = {'--x', '--x'}; % Cross for before

% Define line labels for the plot
line_labels = {'User 1', 'User 2'};

% Create a new figure
figure;
linewidth = 2;
savefig = 1;
fontname = 'Times New Roman';
fontsize = 20;
legendsize = 14;
markersize = 10;
set(gcf, 'position', [0 0 800 600]);
root_save = 'C:\Program Files\UFRJ\TCC\images\';

% Plot BER vs SNR for each user and scenario
for i = 1:length(rows_to_plot)
    semilogy(SNR, ber_after_selected(i, :), markers_after{i}, 'Color', colors(i, :), 'LineWidth', linewidth, 'MarkerSize', markersize, 'DisplayName', sprintf('Project Model - %s', line_labels{i}));
    hold on;
    semilogy(SNR, ber_before_selected(i, :), markers_before{i}, 'Color', colors(i, :), 'LineWidth', linewidth, 'MarkerSize', markersize, 'DisplayName', sprintf('Baseline Model - %s', line_labels{i}));
end

% Customize the plot
xlabel('SNR (in dB)', 'FontName', fontname, 'FontSize', fontsize);
ylabel('BER', 'FontName', fontname, 'FontSize', fontsize);
xlim([-10 35]);
legend('show', 'FontName', fontname, 'FontSize', legendsize, 'Location', 'southwest');
legend box off;

set(gca, 'FontName', fontname, 'FontSize', fontsize);

% Save the figure if savefig is set to 1
if savefig == 1
    saveas(gcf, [root_save 'Comparing_models_BER_vs_SNR'], 'fig');
    saveas(gcf, [root_save 'Comparing_models_BER_vs_SNR'], 'png');
    saveas(gcf, [root_save 'Comparing_models_BER_vs_SNR'], 'epsc2');
end

