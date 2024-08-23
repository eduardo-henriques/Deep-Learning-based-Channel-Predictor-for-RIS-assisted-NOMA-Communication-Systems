% Parameters
P_total = 1;        % Total power (normalized)
B_total = 1;        % Total bandwidth (normalized in MHz)
N0 = 1;             % Noise power spectral density (normalized)
L = 20;             % Number of RIS elements (not used in the current script)

% Simple channel: h1 is equal to h2 for this scenario
h1 = 1;
h2 = 1;

% NOMA Power Allocation
p1 = 0:0.01:1;        % Power allocated to UE1 for NOMA
p2 = 1 - p1;          % Power allocated to UE2 for NOMA
R_NOMA_UE1 = zeros(1, length(p1));
R_NOMA_UE2 = zeros(1, length(p1));

for i = 1:length(p1)
    % Calculate NOMA Throughput with SIC for the stronger user (UE1)
    R_NOMA_UE1(i) = B_total * log2(1 + (p1(i) * P_total * norm(h1)^2) / N0);
    R_NOMA_UE2(i) = B_total * log2(1 + (p2(i) * P_total * norm(h2)^2) / ...
                                   (N0 + p1(i) * P_total * norm(h2)^2));
end

% OMA - Equal Bandwidth Partitioning and Power Allocation
alpha = 0:0.01:1;       % Bandwidth partitioning factor for OMA
R_OMA_UE1 = zeros(1, length(alpha));
R_OMA_UE2 = zeros(1, length(alpha));

for i = 1:length(alpha)
    % Calculate OMA Throughput
    R_OMA_UE1(i) = B_total * alpha(i) * log2(1 + (P_total * norm(h1)^2) / N0);
    R_OMA_UE2(i) = B_total * (1 - alpha(i)) * log2(1 + (P_total * norm(h2)^2) / N0);
end

% Plot Results
savefig = 1;           % Flag to save figures
% root_save = 'C:\Your\Desired\Path\';  % Example path (You might want to write it yourself). Root path to save figures
linewidth = 3;         % Line width for plots
fontname = 'Times New Roman';
fontsize = 20;

legend_alg = {'NOMA', 'OMA'};
location_2 = 'northeast';

colours = [0.0000 0.4470 0.7410;  % Colors for the plots
           0.8500 0.3250 0.0980];

figure;
set(gcf, 'Position', [0 0 800 600]);

plot(R_NOMA_UE1, R_NOMA_UE2, '-', 'Color', colours(1,:), 'LineWidth', linewidth);
hold on;
plot(R_OMA_UE1, R_OMA_UE2, '--', 'Color', colours(2,:), 'LineWidth', linewidth);

xlabel('Throughput of UE_{1} (Mbps)', 'FontName', fontname, 'FontSize', fontsize);
ylabel('Throughput of UE_{2} (Mbps)', 'FontName', fontname, 'FontSize', fontsize);

legend(legend_alg, 'FontName', fontname, 'FontSize', fontsize, 'Location', location_2);
legend boxoff;
set(gca, 'FontName', fontname, 'FontSize', fontsize);

if savefig == 1
    saveas(gcf, [root_save 'noma_oma_throughput_comparison'], 'fig');
    saveas(gcf, [root_save 'noma_oma_throughput_comparison'], 'png');
    saveas(gcf, [root_save 'noma_oma_throughput_comparison'], 'epsc2');
end