clear;
close all;
clc;

% Fix the random seed for reproducibility
seed = 42; % Choose any integer seed value you prefer
rng(seed);

%% Initializing Parameters
K = 20000; % Number of time steps
M = 2; % Number of users
v = 0.01; % Speed of users
L = 20; % Number of RIS elements
% Define a range of SNR values (in dB)
SNR_dB = -10:5:40; % From -10 dB to 30 dB
SNR_linear = 10.^(SNR_dB/10); % Convert dB to linear scale
d0 = 20;
d_sr = 150; 
d_s1 = 30;
d_s2 = 40;
alpha_m = 1;
alpha_r1 = 2.2;
alpha_r2 = 2.2;
alpha_sr = 2.2;
N_t = 1;
N_r = 1;
beta = [0.3, 0.7]; % Power Allocation considering user 1 as the strong user
N_SNR = length(SNR_dB);

%% Generate Dataset
% Constants for Modulation
B     = 4; % Number of bits in a symbol (adjusted for simplicity)
Q     = 2.^B; % Constellation size
N_B   = length(B); % Length of B
m_bit = cell(N_B,M,K,N_SNR);

% Initialization
Y = zeros(2*N_t*N_r, M, K, N_SNR, N_B); % Signals at time step t Y = {Y(1),...,Y(K)} and Y(1) = {y1(1), y2(1)} for M = 2
G = zeros(2*N_t*N_r, M, K, N_SNR, N_B); % Cascaded channel gain in K time steps sequence G = {G(1),...,G(K)} and G(1) = {G1(1), G2(1)} for M = 2
PL = zeros(1, M, K);
% Initialize BER storage for both users under NOMA
BER_1_NOMA = zeros(N_SNR, N_B);
%BER_2_NOMA = zeros(N_SNR, N_B);

% Define beta values
% beta = zeros(1, 2, N_SNR);
% beta(:,:,1:7) = repmat([0.3 0.7], 1, 1, 7);
% beta(:,:,8:10) = repmat([0.4 0.6], 1, 1, 3);
% beta(:,:,11:12) = repmat([0.45 0.55], 1, 1, 2);

d_sj = zeros(K,1);
%SNR_current = SNR_linear(i);
%N0_current = 1/SNR_current; % Assuming unit energy for simplicity

for snr_ind=1:N_SNR
    for t=1:K
        for n=1:N_B
            g_t = zeros(L,M);
            y_t = [];
            G_t = []; 
            h_RIS = zeros(N_r,M);
            w_t = zeros(N_r,M);

            m_bit{n,1,t, snr_ind} = randi([0 1],B(n), 1);
            m_bit{n,2,t, snr_ind} = randi([0 1],B(n), 1);
            m_int_1 = m_bit{n,1,t,snr_ind};
            m_int_2 = m_bit{n,2,t,snr_ind};
            m_int = [m_int_1, m_int_2];
            m_int = bit2int(m_int,B(n));
            m_mod = qammod(m_int,Q(n));                                    % M-QAM modulation
            x = m_mod;
            
            Px = vecnorm(x, 2, 1).^2/size(x,1);                                  % Average transmission power
        
            p1_NOMA = beta(1);
            p2_NOMA = beta(2);
    
            % Channel realization
            h_nLoS = sqrt(0.5)*(randn(M, N_t) + 1j*randn(M, N_t));
            h_LoS = sqrt(0.1*0.5)+1*h_nLoS;
            g_0 = (randn(L, N_t) + 1i*randn(L,N_t))*sqrt(0.1*0.5)+1; % Channel between BS and RIS
            %g = sqrt(0.5)*(randn(M, L) + 1j*randn(M, L))/sqrt(2); % Channel between RIS and each user
        
            %x1_t = (randn(N_t,1) + 1i*randn(N_t,1))*sqrt(0.1*0.5)+1; % Random complex number for x1(t)
            %x2_t = (randn(N_t,1) + 1i*randn(N_t,1))*sqrt(0.1*0.5)+1; % Random complex number for x2(t)
    
            % Superposed NOMA signal
            x_NOMA = sqrt(p1_NOMA/Px(1))*x(:,1) + sqrt(p2_NOMA/Px(2))*x(:,2);
%             x_NOMA(:,1) = sqrt(p1_NOMA/Px(1))*x(:,1);
%             x_NOMA(:,2) = sqrt(p2_NOMA/Px(2))*x(:,2);
%             Px_NOMA(:,1) = norm(x_NOMA(:,1))^2/length(x_NOMA(:,1));
%             Px_NOMA(:,2) = norm(x_NOMA(:,2))^2/length(x_NOMA(:,2));

            Px_NOMA = norm(x_NOMA)^2/length(x_NOMA);

            %x_NOMA = sqrt(1/Px_NOMA)*x_NOMA;
            % Do superposition coding
            %x_t = sqrt(beta(1))*x1_t + sqrt(beta(2))*x2_t;
            %x_t = (randn(N_t,1) + 1i*randn(N_t,1))*sqrt(0.5)+1;
        
            for j = 1:M
            
                % Initializing
                y_jt = zeros(1,2*N_t*N_r);
                G_jt = zeros(1,2*N_t*N_r);
        
                % Generate g_j(t), theta_m(t)
                if j == 1
                    if t == 1
                        d_sj(t) = d_s1;
                    else
                        d_sj(t) = d_sj(t-1) + v;
                    end
                    alpha_rj = alpha_r1;
                    g_t(:,j) = (randn(L, N_r) + 1i*randn(L,N_r))*sqrt(0.5)+4; % Rician fading for each channel between RIS and UE1
                else
                    if t == 1
                        d_sj(t) = d_s2;
                    else
                        d_sj(t) = d_sj(t-1) + v;
                    end
                    alpha_rj = alpha_r2;
                    g_t(:,j) = (randn(L, N_r) + 1i*randn(L,N_r))*sqrt(0.5)+3; % Rician fading for each channel between RIS and UE2
                end
                
                % Direct channel with RIS
                theta = linspace(0.01*pi, 0.02*pi, L); % Phase shift range for L RIS elements
                phi = abs(alpha_m).* exp(1i*theta);
                Theta_matrix = diag(phi);
                h_RIS(j) = g_t(:,j)'*Theta_matrix*g_0;
        
                % Calculate Path Loss PL_j(t)
                PL(:,j,t) = sqrt((d_sr/d0)^(-alpha_sr)) * sqrt((d_sj(t)/d0)^(-alpha_rj));
                %disp(size(g))
                
                w = (randn(N_r,M) + 1i*randn(N_r,M))/sqrt(2); % Gaussian distributed noise with mu = 0 and No = 1
%                 w = (randn(N_r,N_t) + 1i*randn(N_r,N_t))/sqrt(2);

                Pv = vecnorm(w,2,1).^2/size(w,1);                                        % Average noise power
            
                w_t(:,j) = sqrt((Px_NOMA/Pv(:,j))*(10^(-SNR_dB(snr_ind)/10)))*w(:,j);         % True AWGN noise for User 1
                
                %Ps = (10^(SNR_dB(snr_ind)/10));
                % Received signal at user j
                % y_j = PL_jt * h_RIS(j) * r_NOMA + w_t;

                % Calculate received signal y_j(t) at jth UE from equation (2)
                y_jt(1,1) = abs(PL(:,j,t) * h_RIS(j) * x_NOMA + w_t(:,j));
                y_jt(1,2) = angle(PL(:,j,t) * h_RIS(j) * x_NOMA + w_t(:,j));
        
                % Concatenate y_j(t), G_j(t) to yt, Gt
                y_t = [y_t; y_jt];
                % Calculate cascaded channel at jth UE from equation (1)
                G_jt(1,1) = abs(PL(:,j,t) * h_RIS(j));
                G_jt(1,2) = angle(PL(:,j,t) * h_RIS(j));

                G_t = [G_t; G_jt];
            end
            % Concatenate yt, Gt to Y, G
            Y(:,:,t, snr_ind, n) = y_t.';
            % Y(1) = Y(1,:,1)
            G(:,:,t, snr_ind, n) = G_t.';
        end
    end
end
save('dataset_3', "Y","G", "SNR_dB", "Q")
save('transmitted_bit_3', "m_bit", "PL")