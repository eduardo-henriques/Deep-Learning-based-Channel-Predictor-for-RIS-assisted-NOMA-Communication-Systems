%% SISO RIS-Aided NOMA System

Parameters
% Define a range of SNR values (in dB)
SNR_dB = 20; % SNR range in dB
K = 13; % Rician K-factor
numUsers = 2; % Number of NOMA users
MONTE_CARLO = 10000; % Number of Monte Carlo iterations for averaging
P = 1; % Total transmit power
SNR_linear = 10.^(SNR_dB/10);
beta = [0.3, 0.7]; % Power Allocation considering user 1 as the strong user
r = 1; % Amplitude reflection coefficient
L = 20; % Number of RIS elements
W = 1e6;

% Pre-allocate matrices for storing results
avgThroughput_NOMA_LOS = zeros(MONTE_CARLO, length(L));
avgThroughput_NOMA_NLOS = zeros(MONTE_CARLO, length(L));
avgThroughput_NOMA_NLOS_RIS = zeros(MONTE_CARLO, length(L));
avgThroughput_NOMA_LOS_RIS = zeros(MONTE_CARLO, length(L));

for monte_carlo_idx = 1:MONTE_CARLO
    % Loop over SNR values
    for i = 1:length(L)
        SNR_current = SNR_linear;
        N0_current = 1/SNR_current; % Assuming unit energy for simplicity
    
        % Channel realization
        h_nLoS = (randn(numUsers, 1) + 1j*randn(numUsers, 1))/sqrt(2);
        h_LoS = sqrt(K/(K+1)) + sqrt(1/(K+1))*h_nLoS;
        g_0 = sqrt(0.5)*(randn(numUsers, L(i)) + 1j*randn(numUsers, L(i))); % Channel between BS and RIS.
        g = sqrt(0.5)*(randn(numUsers, L(i)) + 1j*randn(numUsers, L(i))); % Channel between RIS and each user

        % Direct channel with RIS
        theta = linspace(0.01*pi, 0.02*pi, L(i)); % Phase shift range for L RIS elements
        phi = abs(r).* exp(1i*theta);
        Theta_matrix = diag(phi);
        h_RIS = h_nLoS + sum(g_0*Theta_matrix*g',2); 
        h_RIS_LOS = h_LoS + sum(g_0*Theta_matrix*g',2);

        % SNR calculations
        % Calculate SNRs for LoS channels
        SNR_1_NOMA_LOS = beta(1) * abs(h_LoS(1)).^2 / N0_current;
        SNR_2_NOMA_LOS = beta(2) * abs(h_LoS(2)).^2 ./ (beta(1) * abs(h_LoS(2)).^2 + N0_current);

        % Calculate SNRs for nLoS channels
        SNR_1_NOMA_NLOS = beta(1) * abs(h_nLoS(1)).^2 / N0_current;
        SNR_2_NOMA_NLOS = beta(2) * abs(h_nLoS(2)).^2 ./ (beta(1) * abs(h_nLoS(2)).^2 + N0_current);

        % Calculate SNRs for nLoS channels with RIS
        SNR_1_NOMA_NLOS_RIS = beta(1) * abs(h_RIS(1)).^2 / N0_current;
        SNR_2_NOMA_NLOS_RIS = beta(2) * abs(h_RIS(2)).^2 ./ (beta(1) * abs(h_RIS(2)).^2 + N0_current);

        % Calculate SNRs for LoS channels with RIS
        SNR_1_NOMA_LOS_RIS = beta(1) * abs(h_RIS_LOS(1)).^2 / N0_current;
        SNR_2_NOMA_LOS_RIS = beta(2) * abs(h_RIS_LOS(2)).^2 ./ (beta(1) * abs(h_RIS_LOS(2)).^2 + N0_current);

        % Calculate capacity
        R_1_NOMA_LOS = W*log2(1 + SNR_1_NOMA_LOS);
        R_2_NOMA_LOS = W*log2(1 + SNR_2_NOMA_LOS);
        R_1_NOMA_NLOS = W*log2(1 + SNR_1_NOMA_NLOS);
        R_2_NOMA_NLOS = W*log2(1 + SNR_2_NOMA_NLOS);
        R_1_NOMA_NLOS_RIS = W*log2(1 + SNR_1_NOMA_NLOS_RIS);
        R_2_NOMA_NLOS_RIS = W*log2(1 + SNR_2_NOMA_NLOS_RIS);
        R_1_NOMA_LOS_RIS = W*log2(1 + SNR_1_NOMA_LOS_RIS);
        R_2_NOMA_LOS_RIS = W*log2(1 + SNR_2_NOMA_LOS_RIS);

        avgThroughput_NOMA_LOS(monte_carlo_idx,i) = (R_1_NOMA_LOS + R_2_NOMA_LOS) / 2;
        avgThroughput_NOMA_NLOS(monte_carlo_idx,i) = (R_1_NOMA_NLOS + R_2_NOMA_NLOS) / 2;
        avgThroughput_NOMA_NLOS_RIS(monte_carlo_idx,i) = (R_1_NOMA_NLOS_RIS + R_2_NOMA_NLOS_RIS) / 2;
        avgThroughput_NOMA_LOS_RIS(monte_carlo_idx,i) = (R_1_NOMA_LOS_RIS + R_2_NOMA_LOS_RIS) / 2;

    end
end

avgThroughput_NOMA_2_LOS = mean(avgThroughput_NOMA_LOS,1);
avgThroughput_NOMA_2_NLOS = mean(avgThroughput_NOMA_NLOS,1);
avgThroughput_NOMA_2_NLOS_RIS = mean(avgThroughput_NOMA_NLOS_RIS,1);
avgThroughput_NOMA_2_LOS_RIS = mean(avgThroughput_NOMA_LOS_RIS,1);

% Data for each scenario
throughputData = {avgThroughput_NOMA_LOS, avgThroughput_NOMA_NLOS, avgThroughput_NOMA_NLOS_RIS};

%% Plot simulation results
% savefig = 0;
% root_save = 'C:\Program Files\UFRJ\TCC\images\';
% linewidth  = 3;
% markersize = 10;
% fontname   = 'Times New Roman';
% fontsize   = 20;
% 
% %legend_alg = {'LoS NOMA', 'nLoS NOMA', 'nLoS RIS-NOMA', 'LoS RIS-NOMA'};
% legend_alg = {'nLoS RIS-NOMA', 'LoS RIS-NOMA'};
% 
% location_1 = 'northwest';
% location_2 = 'northeast';
% location_3 = 'southwest';
% location_4 = 'southeast';
% 
% colours = [0.0000 0.4470 0.7410;
%            0.8500 0.3250 0.0980;
%            0.9290 0.6940 0.1250;
%            0.4940 0.1840 0.5560;
%            0.4660 0.6740 0.1880;
%            0.3010 0.7450 0.9330;
%            0.6350 0.0780 0.1840;
%            0.0000 0.0000 0.0000];
% 
% % Plot ergodic Capacity
% fig1=figure;
% 
% set(gcf,'position',[100 100 800 600]);
% %plot(L, avgThroughput_NOMA_2_LOS*1e-6, 'r-', 'color',colours(1,:),'linewidth',linewidth);
% 
% %plot(L, avgThroughput_NOMA_2_NLOS*1e-6, 'b-','color',colours(2,:),'linewidth',linewidth);
% plot(L, avgThroughput_NOMA_2_NLOS_RIS*1e-6, 'y-','color',colours(3,:),'linewidth',linewidth);
% hold on;
% plot(SNR_dB, avgThroughput_NOMA_2_LOS_RIS*1e-6, 'y-','color',colours(4,:),'linewidth',linewidth);
% 
% xlabel('Number of reflecting elements','fontname',fontname,'fontsize',fontsize);
% ylabel('Average Throughput (in Mbps)','fontname',fontname,'fontsize',fontsize);
% 
% legend(legend_alg,'fontname',fontname,'fontsize',fontsize,'location',location_1);
% legend boxoff;
% set(gca,'fontname',fontname,'fontsize',fontsize);
% 
% if savefig == 1
%     filename = fullfile(root_save, 'RIS_NOMA_vs_n_ris_elements');
%     saveas(fig1, filename, 'fig');
%     saveas(fig1, filename, 'png');
%     saveas(fig1, filename, 'epsc2');
% end
% hold off;


% PLot Outage Probability
% Create the figure for plotting CDFs
fig2 = figure;
set(gcf, 'position', [900, 100, 800, 600]);
hold on;

% Colors corresponding to each scenario
legends = {'LoS NOMA', 'nLoS NOMA', 'nLoS RIS-NOMA'};

% Loop through each scenario and plot the CDF
for i = 1:length(throughputData)
    [f, x] = ecdf(reshape(throughputData{i}, [], 1)); % Compute empirical CDF
    
    % Add a point at the end to ensure the CDF continues at 1
    x = [x; max(x)*1.21]; % Extend x-axis range
    f = [f; 1]; % Add a point at the end with CDF 1
    stairs(x*1e-6, f, 'Color', colours(i,:), 'DisplayName', legends{i}, 'linewidth', linewidth);
end

% Customize the plot
xlabel('Average Throughput (in Mbps)', 'FontName', fontname, 'FontSize', fontsize);
ylabel('Cumulative Distribution Function', 'FontName', fontname, 'FontSize', fontsize);
legend('show', 'FontName', fontname, 'FontSize', fontsize, 'Location', 'southeast');
legend boxoff;
box on;
grid off;
set(gca, 'FontName', fontname, 'FontSize', fontsize);

% Check if the figures should be saved
if savefig == 1
    filename = fullfile(root_save, 'noma_Throughput_CDF_comparison');
    saveas(fig2, filename, 'fig');
    saveas(fig2, filename, 'png');
    saveas(fig2, filename, 'epsc2');
end


