% Parameters
P_total = 1; % Total power for simplicity
B_total = 1; % Total bandwidth in MHz normalized
N0 = 1; % Noise power spectral density normalized
L = 20; % Number of RIS elements

% Simple channel. Channel 1, h1, 20 times greater than h2
h1 = 1;
h2 = 1;

% NOMA Power Allocation
p1 = 0:0.01:1; % Power allocated to UE1 for NOMA
p2 = 1-p1; % Power allocated to UE2 for NOMA
R_NOMA_UE1 = zeros(1, length(p1));
R_NOMA_UE2 = zeros(1, length(p1));

for i = 1:length(p1)
    % Calculate NOMA Rates with SIC for the stronger user (UE1)
    R_NOMA_UE1(i) = B_total * log2(1 + (p1(i) * P_total * norm(h1).^2) / N0);
    R_NOMA_UE2(i) = B_total * log2(1 + (p2(i) * P_total * norm(h2).^2) / (N0 + p1(i) * P_total * norm(h2).^2));
end

% OMA - Equal Bandwidth Partitioning and Power Allocation
alpha = 0:0.01:1; % Bandwidth partitioning factor for OMA
R_OMA_UE1 = zeros(1, length(alpha));
R_OMA_UE2 = zeros(1, length(alpha));

for i = 1:length(alpha)
    % Calculate OMA Rates
    R_OMA_UE1(i) = B_total * alpha(i) * log2(1 + (P_total * norm(h1).^2) / N0);
    R_OMA_UE2(i) = B_total * (1-alpha(i)) * log2(1 + (P_total * norm(h2).^2) / N0);
end

% Plot simulation results
savefig = 1;
root_save = 'C:\Program Files\UFRJ\TCC\images\';
linewidth  = 3;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

legend_alg = {'NOMA', 'OMA'};
location_2 = 'northeast';

colours = [0.0000 0.4470 0.7410;
           0.8500 0.3250 0.0980];

figure;

set(gcf,'position',[0 0 800 600]);

plot(R_NOMA_UE1, R_NOMA_UE2, '-','color',colours(1,:),'linewidth',linewidth);
hold on;
plot(R_OMA_UE1, R_OMA_UE2, '--','color',colours(2,:),'linewidth',linewidth);

xlabel('Throughput of UE_{1} (in Mbps)','fontname',fontname,'fontsize',fontsize);
ylabel('Throughput of UE_{2} (in Mbps)','fontname',fontname,'fontsize',fontsize);

legend(legend_alg,'fontname',fontname,'fontsize',fontsize,'location',location_2);
legend boxoff;
set(gca,'fontname',fontname,'fontsize',fontsize);

if savefig == 1
    saveas(gcf,[root_save 'noma_oma_throughput_comparison_same_channel'], 'fig');
    saveas(gcf,[root_save 'noma_oma_throughput_comparison_same_channel'], 'png');
    saveas(gcf,[root_save 'noma_oma_throughput_comparison_same_channel'], 'epsc2');
end

