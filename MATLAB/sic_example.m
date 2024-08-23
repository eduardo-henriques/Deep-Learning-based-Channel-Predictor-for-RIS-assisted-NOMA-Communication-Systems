Cleaning Workspace
clear;
close all;
clc;
%% Plotting Parameters
linewidth  = 2;
fontname   = 'Times New Roman';
fontsize   = 20;
markersize = 10;

%% Data Transmission Parameters
% Parameters related to the data transmission, e.g., number of bits, constellation size, number of subcarriers, and number of transmitted blocks.
B     = 2;                           % Number of bits in a symbol
Q     = 2.^B;                        % Constellation size
N_BLK = 20000;                       % Number of transmitted data blocks (number of time steps in your case)
M = 2;                               % Number of users

% Environment Parameters
% Parameters related to the environment noise.
SNR   = -10:5:40;                     % Signal-to-noise ratio in dB
N_SNR = length(SNR);                 % Length of SNR

%% Communication System
MONTE_CARLO = 100;                   % One hundred rounds to get the point's average

% Initialize matrices to store the results
y   = zeros(N_BLK,M,N_SNR,MONTE_CARLO);
BER = zeros(M,N_SNR,MONTE_CARLO);



Transmitter
m_bit = randi([0 1],B*N_BLK,M);                         % Random bit generation: B*N_BLK bits for each user
m_mod = qammod(m_bit,Q,'InputType','bit');              % Q-QAM modulation
x = m_mod;

Px = vecnorm(x).^2/size(x,1);                           % Average transmission power

P1 = 0.3;                                               % Power allocated to User 1 (stronger channel)
P2 = 1 - P1;                                            % Power allocated to User 2 (weaker channel)

x_NOMA = sqrt(P1/Px(1))*x(:,1) + sqrt(P2/Px(2))*x(:,2); % Superposed signal

Px_NOMA = norm(x_NOMA)^2/length(x_NOMA);


Receiver
for monte_carlo_idx = 1:MONTE_CARLO
    % Generate channel gains
    h1 = (randn(N_BLK,1) + 1i*randn(N_BLK,1))/sqrt(2) + 4;
    h2 = (randn(N_BLK,1) + 1i*randn(N_BLK,1))/sqrt(2) + 3;
    e1 = (randn(N_BLK,1) + 1i*randn(N_BLK,1))*sqrt(0);
    e2 = (randn(N_BLK,1) + 1i*randn(N_BLK,1))*sqrt(0);
    for snr_ind = 1:N_SNR
        v  = (randn(length(x_NOMA),M) + 1i*randn(length(x_NOMA),M))/sqrt(2); % AWGN noise
        Pv = vecnorm(v).^2/size(v,1);                                        % Average noise power

        v1 = sqrt(((Px_NOMA/Pv(1)))*(10^(-SNR(snr_ind)/10)))*v(:,1);         % True AWGN noise for User 1
        v2 = sqrt(((Px_NOMA/Pv(2)))*(10^(-SNR(snr_ind)/10)))*v(:,2);         % True AWGN noise for User 2

        y(:,1,snr_ind,monte_carlo_idx) = h1.*x_NOMA + v1;                     % Received signal by User 1
        y(:,2,snr_ind,monte_carlo_idx) = h2.*x_NOMA + v2;                     % Received signal by User 2

        % Assuming h_hat = h
        h1_hat = h1 + e1;
        h2_hat = h2 + e2;

        % Decoding for User 2
        x_2 = y(:,1,snr_ind,monte_carlo_idx)./(h1_hat*sqrt(P2/Px(2)));
        e3 = (randn(N_BLK,1) + 1i*randn(N_BLK,1))*sqrt(0);
        x_2_hat = x_2 + e3;

        m_bit_2_hat = qamdemod(x_2_hat,Q,'OutputType','bit');
        m_mod_2_hat = qammod(m_bit_2_hat,Q,'InputType','bit');

        % Decoding for User 1
        x_1_hat = (y(:,1,snr_ind,monte_carlo_idx) - h1_hat.*sqrt(P2/Px(1)).*m_mod_2_hat)./(h1_hat*sqrt(P1/Px(1)));
        m_bit_1_hat = qamdemod(x_1_hat,Q,'OutputType','bit');

        % Direct decode for User 2
        x_2_hat_d = y(:,2,snr_ind,monte_carlo_idx)./(h2_hat*sqrt(P2/Px(2)));
        m_bit_2_hat_d = qamdemod(x_2_hat_d,Q,'OutputType','bit');

        % Calculate BER for each user
        BER(1,snr_ind,monte_carlo_idx) = sum(sum(abs(m_bit(:,1)-m_bit_1_hat(:))))/length(m_bit_1_hat(:));
        BER(2,snr_ind,monte_carlo_idx) = sum(sum(abs(m_bit(:,2)-m_bit_2_hat_d(:))))/length(m_bit_2_hat_d(:));
    end
end

% Average BER over Monte Carlo runs
BER_avg = mean(BER, 3);

% Plotting
figure;
savefig = 1;
set(gcf, 'position', [0 0 800 600]);
root_save = 'C:\Program Files\UFRJ\TCC\images\';

semilogy(SNR,BER_avg(1,:),'linewidth',linewidth);
hold on;
semilogy(SNR,BER_avg(2,:),'linewidth',linewidth);

xlabel('SNR (in dB)','fontname',fontname,'fontsize',fontsize);
ylabel('BER','fontname',fontname,'fontsize',fontsize);
xlim([-10 35]);

legend({'User 1','User 2'},'fontname',fontname,'fontsize',fontsize,'location','southwest')
legend box off

set(gca,'fontname',fontname,'fontsize',fontsize);
if savefig == 1
    saveas(gcf, [root_save 'NOMA_SIC_example'], 'fig');
    saveas(gcf, [root_save 'NOMA_SIC_example'], 'png');
    saveas(gcf, [root_save 'NOMA_SIC_example'], 'epsc2');
end