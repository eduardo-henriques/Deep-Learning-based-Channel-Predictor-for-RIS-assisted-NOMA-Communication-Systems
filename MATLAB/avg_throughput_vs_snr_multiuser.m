% Define a range of SNR values (in dB)
SNR_dB = -10:2:30; % From -10 dB to 30 dB
SNR_linear = 10.^(SNR_dB/10); % Convert dB to linear scale

% System parameters
W = 1e6; % Bandwidth in Hz
alpha = 1/3; % Bandwidth partitioning factor for OMA
P = 1; 

% Adjust power allocation
p1 = 0.1*P; % Power allocated to UE1 (User 1)
p2 = 0.2*P; % Power allocated to UE2 (User 2)
p3 = 0.35*P; % Power allocated to UE3 (User 3)
p4 = 0.35*P; % Power allocated to UE4 (User 4)

MONTE_CARLO = 100; % One hundred rounds to get the point's average

% Initialize arrays to store average throughputs
avgThroughput_NOMA = zeros(MONTE_CARLO, length(SNR_dB));
avgThroughput_OMA = zeros(MONTE_CARLO, length(SNR_dB));

for monte_carlo_idx = 1:MONTE_CARLO
    % Rayleigh fading channels for UE1 and UE2
    h1 = sqrt(80)*(randn + 1i*randn)/sqrt(2); % Channel gain for UE1
    h2 = sqrt(20)*(randn + 1i*randn)/sqrt(2); % Channel gain for UE2
    h3 = (randn + 1i*randn)/sqrt(2); % Channel gain for UE3
    h4 = (randn + 1i*randn)/sqrt(2); % Channel gain for UE4

    for i = 1:length(SNR_dB)
        % Recalculate N0 based on current SNR value
        SNR_current = SNR_linear(i);
        N0_current = 1/SNR_current; % Assuming unit energy for simplicity

        % Calculate SNRs
        SNR_1_NOMA = p1 * abs(h1).^2 / N0_current;
        SNR_2_NOMA = p2 * abs(h2).^2 ./ (p1 * abs(h2).^2 + N0_current);
        SNR_3_NOMA = p3 * abs(h3).^2 ./ (p1 * abs(h3).^2 + p2 * abs(h3).^2+ N0_current);
        SNR_4_NOMA = p4 * abs(h4).^2 ./ (p1 * abs(h4).^2 + p2 * abs(h4).^2 + p3 * abs(h4).^2 + N0_current);
        SNR_1_OMA = (P/3) * abs(h1).^2 / (alpha * N0_current);
        SNR_2_OMA = (P/3) * abs(h2).^2 / ((alpha) * N0_current);
        SNR_3_OMA = (P/3) * abs(h3).^2 / ((alpha) * N0_current);
        SNR_4_OMA = P * abs(h4).^2 / ((alpha) * N0_current);

        R_1_NOMA = W * log2(1 + SNR_1_NOMA);
        R_2_NOMA = W * log2(1 + SNR_2_NOMA);
        R_3_NOMA = W * log2(1 + SNR_3_NOMA);
        R_4_NOMA = W * log2(1 + SNR_4_NOMA);
        R_1_OMA = W * alpha * log2(1 + SNR_1_OMA);
        R_2_OMA = W * (alpha) * log2(1 + SNR_2_OMA);
        R_3_OMA = W * (alpha) * log2(1 + SNR_3_OMA);
        R_4_OMA = W * (alpha) * log2(1 + SNR_4_OMA);

        % Calculate average throughputs for NOMA and OMA
        avgThroughput_NOMA(monte_carlo_idx,i) = (R_1_NOMA + R_2_NOMA + R_3_NOMA + R_4_NOMA) / 4;
        avgThroughput_OMA(monte_carlo_idx,i) = (R_1_OMA + R_2_OMA + R_3_OMA + R_4_NOMA) / 4;
    end
end

avgThroughput_NOMA_2 = mean(avgThroughput_NOMA,1);
avgThroughput_OMA_2 = mean(avgThroughput_OMA,1);

% % Plotting the results
% figure;
% plot(SNR_dB, avgThroughput_NOMA_2, 'r-', 'LineWidth', 2);
% hold on;
% plot(SNR_dB, avgThroughput_OMA_2, 'b-', 'LineWidth', 2);
% plot(SNR_dB, avgThroughput_NOMA_2, 'r--', 'LineWidth', 2);
% plot(SNR_dB, avgThroughput_OMA_2, 'b--', 'LineWidth', 2);
% 
% xlabel('SNR (dB)');
% ylabel('Average Throughput (bps)');
% title('Average Throughput vs. SNR for NOMA and OMA with Rayleigh Fading');
% legend('NOMA - p_1 = 0.4', 'OMA - \alpha = 0.5','NOMA - p_1 = 0.7', 'OMA - \alpha = 0.6');
% grid on;

% Plot simulation results
savefig = 1;
root_save = 'C:\Program Files\UFRJ\TCC\images\';
linewidth  = 3;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

legend_alg = {'OMA', 'NOMA'};

location_1 = 'northwest';
location_2 = 'northeast';
location_3 = 'southwest';
location_4 = 'southeast';

colours = [0.0000 0.4470 0.7410;
           0.8500 0.3250 0.0980;
           0.9290 0.6940 0.1250;
           0.4940 0.1840 0.5560;
           0.4660 0.6740 0.1880;
           0.3010 0.7450 0.9330;
           0.6350 0.0780 0.1840;
           0.0000 0.0000 0.0000];

figure;

set(gcf,'position',[0 0 800 600]);
plot(SNR_dB, avgThroughput_OMA_2, 'r-', 'color',colours(2,:),'linewidth',linewidth);
hold on;
plot(SNR_dB, avgThroughput_NOMA_2, 'b-','color',colours(1,:),'linewidth',linewidth);

xlabel('SNR (in dB)','fontname',fontname,'fontsize',fontsize);
ylabel('Average Throughput (in bps)','fontname',fontname,'fontsize',fontsize);

legend(legend_alg,'fontname',fontname,'fontsize',fontsize,'location',location_4);
legend boxoff;
set(gca,'fontname',fontname,'fontsize',fontsize);

if savefig == 1
    saveas(gcf,[root_save 'noma_oma_average_throughput_vs_snr_comparison_3_users'], 'fig');
    saveas(gcf,[root_save 'noma_oma_average_throughput_vs_snr_comparison_3_users'], 'png');
    saveas(gcf,[root_save 'noma_oma_average_throughput_vs_snr_comparison_3_users'], 'epsc2');
end


