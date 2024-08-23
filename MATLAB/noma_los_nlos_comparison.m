%% Setting the Parameters
SNR_dB = -10:2:30;             % SNR range in dB
SNR_linear = 10.^(SNR_dB/10);  % Convert SNR from dB to linear scale
rng(42);  % Set seed for reproducibility

% System parameters
W = 1e6;        % Bandwidth in Hz
P = 1;          % Total power (normalized)
N0 = 1e-21;     % Noise power spectral density (normalized)
beta = [0.3, 0.7];  % Power allocation for User 1 and User 2, respectively
num_users = 2;

%% Channel Model
% nLoS channel gains (Rayleigh fading)
h_nLoS = (randn(1, num_users) + 1i*randn(1, num_users)) / sqrt(2);

% LoS channel gains (Rician fading with K-factor)
K = 13;  % Rician K-factor
h_LoS = sqrt(K/(K+1)) + sqrt(1/(K+1)) * h_nLoS;

%% Simulation Settings
MONTE_CARLO = 10000;  % Number of Monte Carlo simulations

% Initialize arrays to store average throughputs
avgThroughput_NOMA_LOS = zeros(MONTE_CARLO, length(SNR_dB));
avgThroughput_NOMA_NLOS = zeros(MONTE_CARLO, length(SNR_dB));

for mc_idx = 1:MONTE_CARLO
    for i = 1:length(SNR_dB)
        % Calculate noise power based on current SNR
        SNR_current = SNR_linear(i);
        N0_current = 1 / SNR_current;

        % SNRs for LoS channels
        SNR_1_NOMA_LOS = beta(1) * abs(h_LoS(1)).^2 / N0_current;
        SNR_2_NOMA_LOS = beta(2) * abs(h_LoS(2)).^2 / (beta(1) * abs(h_LoS(2)).^2 + N0_current);

        % SNRs for nLoS channels
        SNR_1_NOMA_NLOS = beta(1) * abs(h_nLoS(1)).^2 / N0_current;
        SNR_2_NOMA_NLOS = beta(2) * abs(h_nLoS(2)).^2 / (beta(1) * abs(h_nLoS(2)).^2 + N0_current);

        % Channel capacity for LoS channels
        R_1_NOMA_LOS = W * log2(1 + SNR_1_NOMA_LOS);
        R_2_NOMA_LOS = W * log2(1 + SNR_2_NOMA_LOS);

        % Channel capacity for nLoS channels
        R_1_NOMA_NLOS = W * log2(1 + SNR_1_NOMA_NLOS);
        R_2_NOMA_NLOS = W * log2(1 + SNR_2_NOMA_NLOS);

        % Calculate average throughputs
        avgThroughput_NOMA_LOS(mc_idx, i) = (R_1_NOMA_LOS + R_2_NOMA_LOS) / 2;
        avgThroughput_NOMA_NLOS(mc_idx, i) = (R_1_NOMA_NLOS + R_2_NOMA_NLOS) / 2;
    end
end

% Compute the mean throughput over Monte Carlo simulations
avgThroughput_NOMA_2_LOS = mean(avgThroughput_NOMA_LOS, 1);
avgThroughput_NOMA_2_NLOS = mean(avgThroughput_NOMA_NLOS, 1);

%% Plotting Results

linewidth  = 3;
fontname   = 'Times New Roman';
fontsize   = 20;

legend_alg = {'NOMA LoS', 'NOMA nLoS'};

figure;
set(gcf, 'Position', [0 0 800 600]);

plot(SNR_dB, avgThroughput_NOMA_2_LOS * 1e-6, '-', 'Color', [0.0000 0.4470 0.7410], 'LineWidth', linewidth);
hold on;
plot(SNR_dB, avgThroughput_NOMA_2_NLOS * 1e-6, '--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', linewidth);

xlabel('SNR (in dB)', 'FontName', fontname, 'FontSize', fontsize);
ylabel('Average Throughput (in Mbps)', 'FontName', fontname, 'FontSize', fontsize);

legend(legend_alg, 'FontName', fontname, 'FontSize', fontsize, 'Location', 'northwest');
legend boxoff;
set(gca, 'FontName', fontname, 'FontSize', fontsize);

% Save figures if required
savefig = 0;
% root_save = 'C:\Your\Desired\Path\';  Example path (You might want to write it yourself). Root path to save figures

if savefig == 1
    saveas(gcf, [root_save 'noma_los_nlos_average_throughput_vs_snr_comparison'], 'fig');
    saveas(gcf, [root_save 'noma_los_nlos_average_throughput_vs_snr_comparison'], 'png');
    saveas(gcf, [root_save 'noma_los_nlos_average_throughput_vs_snr_comparison'], 'epsc2');
end