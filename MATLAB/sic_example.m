% Clean Workspace
clear;
close all;
clc;

%% Plotting Parameters
linewidth  = 2;
fontname   = 'Times New Roman';
fontsize   = 20;
markersize = 10;

%% Data Transmission Parameters
B     = 2;                  % Number of bits per symbol
Q     = 2^B;                % Constellation size (Q-QAM)
N_BLK = 20000;              % Number of transmitted data blocks (time steps)
M     = 2;                  % Number of users

%% Environment Parameters
SNR   = -10:5:40;           % Signal-to-noise ratio in dB
N_SNR = length(SNR);        % Number of SNR points

%% Simulation Settings
MONTE_CARLO = 100;          % Number of Monte Carlo runs

% Initialize result storage
y   = zeros(N_BLK, M, N_SNR, MONTE_CARLO);
BER = zeros(M, N_SNR, MONTE_CARLO);

%% Transmitter
m_bit = randi([0 1], B*N_BLK, M);                % Random bit generation for each user
m_mod = qammod(m_bit, Q, 'InputType', 'bit');    % Q-QAM modulation

% Power allocation
P1 = 0.3;                                        % Power allocated to User 1 (stronger channel)
P2 = 1 - P1;                                     % Power allocated to User 2 (weaker channel)

x_NOMA = sqrt(P1)*m_mod(:,1) + sqrt(P2)*m_mod(:,2); % Superposed NOMA signal

%% Receiver
for mc_idx = 1:MONTE_CARLO
    % Channel gains for each user
    h1 = (randn(N_BLK,1) + 1i*randn(N_BLK,1))/sqrt(2) + 4;
    h2 = (randn(N_BLK,1) + 1i*randn(N_BLK,1))/sqrt(2) + 3;

    for snr_idx = 1:N_SNR
        % AWGN noise
        v = (randn(N_BLK, M) + 1i*randn(N_BLK, M))/sqrt(2); 
        
        % Calculate noise variance based on SNR
        v1 = sqrt(10^(-SNR(snr_idx)/10)) * v(:,1);
        v2 = sqrt(10^(-SNR(snr_idx)/10)) * v(:,2);

        % Received signals
        y(:,1,snr_idx,mc_idx) = h1 .* x_NOMA + v1;
        y(:,2,snr_idx,mc_idx) = h2 .* x_NOMA + v2;

        % SIC for User 1
        x_2_hat = y(:,1,snr_idx,mc_idx) ./ (h1 * sqrt(P2));
        m_bit_2_hat = qamdemod(x_2_hat, Q, 'OutputType', 'bit');
        m_mod_2_hat = qammod(m_bit_2_hat, Q, 'InputType', 'bit');

        x_1_hat = (y(:,1,snr_idx,mc_idx) - h1 .* sqrt(P2) .* m_mod_2_hat) ./ (h1 * sqrt(P1));
        m_bit_1_hat = qamdemod(x_1_hat, Q, 'OutputType', 'bit');

        % Direct decode for User 2
        x_2_hat_d = y(:,2,snr_idx,mc_idx) ./ (h2 * sqrt(P2));
        m_bit_2_hat_d = qamdemod(x_2_hat_d, Q, 'OutputType', 'bit');

        % BER Calculation
        BER(1,snr_idx,mc_idx) = sum(sum(m_bit(:,1) ~= m_bit_1_hat)) / numel(m_bit_1_hat);
        BER(2,snr_idx,mc_idx) = sum(sum(m_bit(:,2) ~= m_bit_2_hat_d)) / numel(m_bit_2_hat_d);
    end
end

% Average BER over Monte Carlo runs
BER_avg = mean(BER, 3);

%% Plotting Results
figure;
set(gcf, 'Position', [0 0 800 600]);

semilogy(SNR, BER_avg(1,:), 'LineWidth', linewidth);
hold on;
semilogy(SNR, BER_avg(2,:), 'LineWidth', linewidth);

xlabel('SNR (in dB)', 'FontName', fontname, 'FontSize', fontsize);
ylabel('BER', 'FontName', fontname, 'FontSize', fontsize);
xlim([-10 35]);

legend({'User 1', 'User 2'}, 'FontName', fontname, 'FontSize', fontsize, 'Location', 'southwest');
legend boxoff;

set(gca, 'FontName', fontname, 'FontSize', fontsize);

% Save figures if required
savefig = 0;
% root_save = 'C:\Your\Desired\Path\';  Example path (You might want to write it yourself). Root path to save figures

if savefig == 1
    saveas(gcf, [root_save 'NOMA_SIC_example'], 'fig');
    saveas(gcf, [root_save 'NOMA_SIC_example'], 'png');
    saveas(gcf, [root_save 'NOMA_SIC_example'], 'epsc2');
end