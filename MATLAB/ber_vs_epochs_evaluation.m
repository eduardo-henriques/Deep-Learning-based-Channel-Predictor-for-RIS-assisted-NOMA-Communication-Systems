clear;
clc;

% Define the epoch values and corresponding file names
epochs = [20, 50, 100];
files = {'ber_epochs_20.mat', 'ber_epochs_50.mat', 'ber_epochs_100.mat'};

% Define the SNR values
SNR_dB = -10:5:40;

% Initialize a cell array to store the BER values
BER_epochs = cell(length(epochs), 1);

% Load each file and extract the BER values
for i = 1:length(epochs)
    load(files{i});
    BER_epochs{i} = BER; % Store the BER matrix for the current epoch
end

% Define the specific SNR value to plot
chosen_snr_value = 30; % SNR of 20 dB
chosen_snr_index = find(SNR_dB == chosen_snr_value);

% Define the BER lines to plot and their labels
ber_lines = [3, 7]; % User 1 and User 2
ber_line_labels = {'User 1', 'User 2'};

% Define colors for each user
user_colors = [0.0000 0.4470 0.7410;  % Blue for User 1
               0.8500 0.3250 0.0980]; % Red for User 2

% Plot settings
savefig = 1;
root_save = 'C:\Program Files\UFRJ\TCC\images\';
linewidth = 2;
markersize = 10;
fontname = 'Times New Roman';
fontsize = 20;

% Create a figure for the plot
figure;
set(gcf, 'position', [0 0 800 600]);
hold on;

% Plot BER vs. epochs for the chosen SNR value and each chosen BER line
for line_idx = 1:length(ber_lines)
    ber_line = ber_lines(line_idx);
    BER_values = arrayfun(@(epoch_idx) BER_epochs{epoch_idx}(ber_line, chosen_snr_index), 1:length(epochs));
    semilogy(epochs, BER_values, '-o', 'Color', user_colors(line_idx, :), 'LineWidth', linewidth, 'MarkerSize', markersize, 'DisplayName', ber_line_labels{line_idx});
end

xlabel('Epochs', 'FontName', fontname, 'FontSize', fontsize);
ylabel('BER', 'FontName', fontname, 'FontSize', fontsize);
grid off;
box on;
legend('show', 'FontName', fontname, 'FontSize', fontsize, 'Location', 'northwest');
legend boxoff;

set(gca, 'FontName', fontname, 'FontSize', fontsize);
hold off;

% Save the figure if required
if savefig == 1
    saveas(gcf, fullfile(root_save, 'BER_vs_Epochs_new_model'), 'fig');
    saveas(gcf, fullfile(root_save, 'BER_vs_Epochs_new_model'), 'png');
    saveas(gcf, fullfile(root_save, 'BER_vs_Epochs_new_model'), 'epsc2');
end

