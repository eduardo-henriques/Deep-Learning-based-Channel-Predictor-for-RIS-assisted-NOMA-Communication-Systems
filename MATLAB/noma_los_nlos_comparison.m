%% Setting the parameters
% Define a range of SNR values (in dB)
SNR_dB = -10:2:30; % From -10 dB to 30 dB
SNR_linear = 10.^(SNR_dB/10); % Convert dB to linear scale

% System parameters
W = 1e6; % Bandwidth in Hz
P = 1; % Total power
N0 = 1e-21; % Noise power spectral density
% Power Allocation Coefficients for NOMA
beta = [0.3, 0.7]; % Power allocated to User 1 and User 2, respectively
num_users = 2;

%% Channel Model
% nLoS channel gains
h_nLoS = (randn(1,2) + 1i*randn(1,2))/sqrt(2); % Example path losses for 2 users at distances 1 and 2 units
% Example Rician fading for LoS, K-factor represents the ratio of direct to scattered components
K = 13; % Example K-factor for LoS condition
h_LoS = sqrt(K/(K+1)) + sqrt(1/(K+1))*h_nLoS;


Model
MONTE_CARLO = 500; % One hundred rounds to get the point's average

% Initialize arrays to store average throughputs
avgThroughput_NOMA_LOS = zeros(MONTE_CARLO, length(SNR_dB));
avgThroughput_NOMA_NLOS = zeros(MONTE_CARLO, length(SNR_dB));

for monte_carlo_idx = 1:MONTE_CARLO
    for i = 1:length(SNR_dB)
        % Recalculate N0 based on current SNR value
        SNR_current = SNR_linear(i);
        N0_current = 1/SNR_current; % Assuming unit energy for simplicity

        % Calculate SNRs for LoS channels
        SNR_1_NOMA_LOS = beta(1) * abs(h_LoS(1)).^2 / N0_current;
        SNR_2_NOMA_LOS = beta(2) * abs(h_LoS(2)).^2 ./ (beta(1) * abs(h_LoS(2)).^2 + N0_current);

        % Calculate SNRs for nLoS channels
        SNR_1_NOMA_NLOS = beta(1) * abs(h_nLoS(1)).^2 / N0_current;
        SNR_2_NOMA_NLOS = beta(2) * abs(h_nLoS(2)).^2 ./ (beta(1) * abs(h_nLoS(2)).^2 + N0_current);

        %Calculate channel capacity for LoS channels
        R_1_NOMA_LOS = W * log2(1 + SNR_1_NOMA_LOS);
        R_2_NOMA_LOS = W * log2(1 + SNR_2_NOMA_LOS);


        R_1_NOMA_NLOS = W * log2(1 + SNR_1_NOMA_NLOS);
        R_2_NOMA_NLOS = W * log2(1 + SNR_2_NOMA_NLOS);

        % Calculate average throughputs for NOMA and OMA
        avgThroughput_NOMA_LOS(monte_carlo_idx,i) = (R_1_NOMA_LOS + R_2_NOMA_LOS) / 2;
        avgThroughput_NOMA_NLOS(monte_carlo_idx,i) = (R_1_NOMA_NLOS + R_2_NOMA_NLOS) / 2;
    end
end

avgThroughput_NOMA_2_LOS = mean(avgThroughput_NOMA_LOS,1);
avgThroughput_NOMA_2_NLOS = mean(avgThroughput_NOMA_NLOS,1);


%% Plot simulation results
savefig = 1;
root_save = 'C:\Program Files\UFRJ\TCC\images\';
linewidth  = 3;
markersize = 10;
fontname   = 'Times New Roman';
fontsize   = 20;

legend_alg = {'NOMA LoS', 'NOMA nLoS'};

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
plot(SNR_dB, avgThroughput_NOMA_2_LOS*1e-6, 'r-', 'color',colours(1,:),'linewidth',linewidth);
hold on;
plot(SNR_dB, avgThroughput_NOMA_2_NLOS*1e-6, 'b-','color',colours(2,:),'linewidth',linewidth);

xlabel('SNR (in dB)','fontname',fontname,'fontsize',fontsize);
ylabel('Average Throughput (in Mbps)','fontname',fontname,'fontsize',fontsize);

legend(legend_alg,'fontname',fontname,'fontsize',fontsize,'location',location_1);
legend boxoff;
set(gca,'fontname',fontname,'fontsize',fontsize);

if savefig == 1
    saveas(gcf,[root_save 'noma_los_nlos_average_throughput_vs_snr_comparison'], 'fig');
    saveas(gcf,[root_save 'noma_los_nlos_average_throughput_vs_snr_comparison'], 'png');
    saveas(gcf,[root_save 'noma_los_nlos_average_throughput_vs_snr_comparison'], 'epsc2');
end
