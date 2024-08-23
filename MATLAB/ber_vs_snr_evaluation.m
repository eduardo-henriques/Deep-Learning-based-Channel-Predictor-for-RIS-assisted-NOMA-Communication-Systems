clc;
clear;

% Load the results
load('./new_model/results_2_users_dl_4QAM_7.mat');
load('transmitted_bit_1.mat');

new_m_bit = m_bit(:,:,21:end,:);
test_m_bit = m_bit(:,:,16005:end,:);
train_m_bit = m_bit(:,:,21:16004,:);


% Function to build a signal from amplitude and phase
build_signal = @(amplitude, phase) abs(amplitude) .* exp(1i * phase);

% Function to calculate Bit Error Rate (BER)
calculate_ber = @(demodulated, original) sum(demodulated ~= original) / numel(original);

num_samples = size(y_pred_all, 1);
test_num_samples = size(y_pred_test, 1);
train_num_samples = size(y_pred_train, 1);
Ps = 10.^(double(SNR_dB)/10)*1e-3;
num_snr = length(SNR_dB);
num_mod = length(Q);
M = size(y_pred_all, 3);  % Number of users
BER = zeros(10, num_snr, num_mod);
P1 = 0.3;
P2 = 1 - P1;

for snr_idx = 1:num_snr
    for mod_idx = 1:num_mod
        ber_1 = zeros(num_samples-1, 1);
        ber_2 = zeros(num_samples-1, 1);
        ber_3 = zeros(test_num_samples-1, 1);
        ber_4 = zeros(test_num_samples-1, 1);
        ber_13 = zeros(test_num_samples-1, 1);
        ber_14 = zeros(test_num_samples-1, 1);
        ber_5 = zeros(train_num_samples-1, 1);
        ber_6 = zeros(train_num_samples-1, 1);
        ber_15 = zeros(train_num_samples-1, 1);
        ber_16 = zeros(train_num_samples-1, 1);

        %% All curve
%         for sample_pos = 2:num_samples
%             received_amplitude_1 = X_all(sample_pos, 20, 1, 1, snr_idx, mod_idx);
%             received_phase_1 = X_all(sample_pos, 20, 1, 2, snr_idx, mod_idx);
%             y_1 = build_signal(received_amplitude_1, received_phase_1);
% 
% %           received_amplitude_2 = X_test_all(sample_pos, 20, 2, 1, snr_idx);
% %           received_phase_2 = X_test_all(sample_pos, 20, 2, 2, snr_idx);
% %           y_2 = build_signal(received_amplitude_2, received_phase_2);
% 
%             pred_channel_amplitude = y_pred_all(sample_pos-1, 1, :, 1, snr_idx, mod_idx);
%             pred_channel_phase = y_pred_all(sample_pos-1, 1, :, 2, snr_idx, mod_idx);
%             h1_hat = build_signal(pred_channel_amplitude(:,:,1,:,:), pred_channel_phase(:,:,1,:,:));
% %           h2_hat = build_signal(pred_channel_amplitude(:,:,2,:,:), pred_channel_phase(:,:,2,:,:));
% 
% 
%             true_channel_amplitude = y_true_all(sample_pos-1, 1, :, 1, snr_idx, mod_idx);
%             true_channel_phase = y_true_all(sample_pos-1, 1, :, 2, snr_idx, mod_idx);
%             h1_true = build_signal(true_channel_amplitude(:,:,1,:,:), true_channel_phase(:,:,1,:,:));
% %           h2_true = build_signal(true_channel_amplitude(:,:,2,:,:), true_channel_phase(:,:,2,:,:));
% 
% %           x_2_hat = y_2./(h2_hat);
% %           b_2_hat = qamdemod(x_2_hat,Q,'OutputType','bit');
%             %m_mod_2_hat = qammod(b_2_hat,double(Q),'InputType','bit');
% 
%             %x_1_hat = (y_1 - h1_hat.*sqrt(P2).*m_mod_2_hat)./(h1_hat*sqrt(P1));
%             x_1_hat = y_1./(h1_hat);
% 
%             x_1_true = y_1./(h1_true);
%             
%             %x_original = estimate_sent_signal(received_signal, true_channel_amplitude, true_channel_phase);
% 
%             b_1_hat = qamdemod(x_1_hat, Q, 'OutputType','bit');
%             b_1_true = qamdemod(x_1_true, Q, 'OutputType','bit');
%             %b_original = qamdemod(s_original, Q);
% 
%             ber_1(sample_pos-1) = calculate_ber(b_1_hat, new_m_bit(mod_idx,1,sample_pos-1, snr_idx));
%             ber_2(sample_pos-1) = calculate_ber(b_1_true, new_m_bit(mod_idx,1,sample_pos-1, snr_idx));
%         end
% 
%         %% Test curve
        for test_num_sample = 2:test_num_samples
            received_test_amplitude_1 = X_test(test_num_sample, 20, 1, 1, snr_idx, mod_idx);
            received_test_phase_1 = X_test(test_num_sample, 20, 1, 2, snr_idx, mod_idx);
            y_1_test = build_signal(received_test_amplitude_1, received_test_phase_1);

            received_test_amplitude_2 = X_test(test_num_sample, 20, 2, 1, snr_idx, mod_idx);
            received_test_phase_2 = X_test(test_num_sample, 20, 2, 2, snr_idx, mod_idx);
            y_2_test = build_signal(received_test_amplitude_2, received_test_phase_2);

            pred_test_channel_amplitude = y_pred_test(test_num_sample-1, 1, :, 1, snr_idx, mod_idx);
            pred_test_channel_phase = y_pred_test(test_num_sample-1, 1, :, 2, snr_idx, mod_idx);
            h1_test_hat = build_signal(pred_test_channel_amplitude(:,:,1,:,:), pred_test_channel_phase(:,:,1,:,:));
            h2_test_hat = build_signal(pred_test_channel_amplitude(:,:,2,:,:), pred_test_channel_phase(:,:,2,:,:));


            true_test_channel_amplitude = y_true_test(test_num_sample-1, 1, :, 1, snr_idx, mod_idx);
            true_test_channel_phase = y_true_test(test_num_sample-1, 1, :, 2, snr_idx, mod_idx);
            h1_test_true = build_signal(true_test_channel_amplitude(:,:,1,:,:), true_test_channel_phase(:,:,1,:,:));
            h2_test_true = build_signal(true_test_channel_amplitude(:,:,2,:,:), true_test_channel_phase(:,:,2,:,:));

            % Decodifica x2 no sinal de y1
            x_2_test_hat = y_1_test./(h1_test_hat*sqrt(P2/2));
            b_2_test_hat = qamdemod(x_2_test_hat,double(Q(mod_idx)),'OutputType','bit');
            m_mod_2_hat = qammod(b_2_test_hat,double(double(Q(mod_idx))),'InputType','bit');

            % Decodifica x2 no sinal de y1
            x_2_test_true = y_1_test./(h1_test_true*sqrt(P2/2));
            b_2_test_true = qamdemod(x_2_test_true, double(Q(mod_idx)), 'OutputType','bit');
            m_mod_2_true = qammod(b_2_test_true,double(double(Q(mod_idx))),'InputType','bit');
            
            % Decodifica x1 
            x_1_test_hat = (y_1_test - h1_test_hat.*sqrt(P2/2).*m_mod_2_hat)./h1_test_hat*sqrt(P1/2);
%             x_1_test_hat = y_1_test./(h1_test_hat);
            x_1_test_true = (y_1_test - h1_test_true.*sqrt(P2/2).*m_mod_2_true)./h1_test_true*sqrt(P1/2);
%             x_1_test_true = y_1_test./(h1_test_true);
            
            % Direct decoding de x2 em y2
            % Decodifica x2 no sinal de y1
            x_2_test_hat_d = y_2_test./(h2_test_hat*sqrt(P2/2));
            b_2_test_hat_d = qamdemod(x_2_test_hat_d,double(Q(mod_idx)),'OutputType','bit');

            % Decodifica x2 no sinal de y1
            x_2_test_true_d = y_2_test./(h2_test_true*sqrt(P2/2));
            b_2_test_true_d = qamdemod(x_2_test_true_d, double(Q(mod_idx)), 'OutputType','bit');
            %x_original = estimate_sent_signal(received_signal, true_channel_amplitude, true_channel_phase);

            b_1_test_hat = qamdemod(x_1_test_hat, double(Q(mod_idx)), 'OutputType','bit');
            b_1_test_true = qamdemod(x_1_test_true, double(Q(mod_idx)), 'OutputType','bit');
            %b_original = qamdemod(s_original, Q);
            target_b_1 = test_m_bit{mod_idx,1,test_num_sample-1, snr_idx};
            target_b_2 = test_m_bit{mod_idx,2,test_num_sample-1, snr_idx};

            ber_3(test_num_sample-1) = calculate_ber(b_1_test_hat, target_b_1);
            ber_4(test_num_sample-1) = calculate_ber(b_1_test_true, target_b_1);
            ber_13(test_num_sample-1) = calculate_ber(b_2_test_hat_d, target_b_2);
            ber_14(test_num_sample-1) = calculate_ber(b_2_test_true_d, target_b_2);
        end

%         %% Train curve
        for train_num_sample = 2:train_num_samples
            received_train_amplitude_1 = X_train(train_num_sample, 20, 1, 1, snr_idx, mod_idx);
            received_train_phase_1 = X_train(train_num_sample, 20, 1, 2, snr_idx, mod_idx);
            y_1_train = build_signal(received_train_amplitude_1, received_train_phase_1);

            received_train_amplitude_2 = X_train(train_num_sample, 20, 2, 1, snr_idx, mod_idx);
            received_train_phase_2 = X_train(train_num_sample, 20, 2, 2, snr_idx, mod_idx);
            y_2_train = build_signal(received_train_amplitude_2, received_train_phase_2);

            pred_train_channel_amplitude = y_pred_train(train_num_sample-1, 1, :, 1, snr_idx, mod_idx);
            pred_train_channel_phase = y_pred_train(train_num_sample-1, 1, :, 2, snr_idx, mod_idx);
            h1_train_hat = build_signal(pred_train_channel_amplitude(:,:,1,:,:), pred_train_channel_phase(:,:,1,:,:));
            h2_train_hat = build_signal(pred_train_channel_amplitude(:,:,2,:,:), pred_train_channel_phase(:,:,2,:,:));


            true_train_channel_amplitude = y_true_train(train_num_sample-1, 1, :, 1, snr_idx, mod_idx);
            true_train_channel_phase = y_true_train(train_num_sample-1, 1, :, 2, snr_idx, mod_idx);
            h1_train_true = build_signal(true_train_channel_amplitude(:,:,1,:,:), true_train_channel_phase(:,:,1,:,:));
            h2_train_true = build_signal(true_train_channel_amplitude(:,:,2,:,:), true_train_channel_phase(:,:,2,:,:));

            % Decodifica x2 em y1
            x_2_train_hat = y_1_train./(h1_train_hat*sqrt(P2/2));
            b_2_train_hat = qamdemod(x_2_train_hat,double(Q(mod_idx)),'OutputType','bit');
            m_mod_2_train_hat = qammod(b_2_train_hat,double(double(Q(mod_idx))),'InputType','bit');

            x_2_train_true = y_1_train./(h1_train_true*sqrt(P2/2));
            b_2_train_true = qamdemod(x_2_train_true, double(Q(mod_idx)), 'OutputType','bit');
            m_mod_2_train_true = qammod(b_2_train_true,double(double(Q(mod_idx))),'InputType','bit');

            x_1_train_hat = (y_1_train - h1_train_hat.*sqrt(P2/2).*m_mod_2_train_hat)./h1_train_hat*sqrt(P1/2);
            %x_1_train_hat = y_1_train./(h1_train_hat);

            %x_1_train_true = y_1_train./(h1_train_true);
            x_1_train_true  = (y_1_train - h1_train_true.*sqrt(P2/2).*m_mod_2_train_true)./h1_train_true*sqrt(P1/2);
            %x_original = estimate_sent_signal(received_signal, true_channel_amplitude, true_channel_phase);

            % Direct decoding de x2 em y2
            x_2_train_hat_d = y_2_train./(h2_train_hat*sqrt(P2/2));
            b_2_train_hat_d = qamdemod(x_2_train_hat_d,double(Q(mod_idx)),'OutputType','bit');

            x_2_train_true_d = y_2_train./(h2_train_true*sqrt(P2/2));
            b_2_train_true_d = qamdemod(x_2_train_true_d, double(Q(mod_idx)), 'OutputType','bit');

            b_1_train_hat = qamdemod(x_1_train_hat, double(Q(mod_idx)), 'OutputType','bit');
            b_1_train_true = qamdemod(x_1_train_true, double(Q(mod_idx)), 'OutputType','bit');
            %b_original = qamdemod(s_original, Q);
            target_b_1_train = train_m_bit{mod_idx,1,train_num_sample-1, snr_idx};
            target_b_2_train = train_m_bit{mod_idx,2,train_num_sample-1, snr_idx};

            ber_5(train_num_sample-1) = calculate_ber(b_1_train_hat, target_b_1_train);
            ber_6(train_num_sample-1) = calculate_ber(b_1_train_true, target_b_1_train);
            ber_15(train_num_sample-1) = calculate_ber(b_2_train_hat_d, target_b_2_train);
            ber_16(train_num_sample-1) = calculate_ber(b_2_train_true_d, target_b_2_train);
        end
        
%         BER(1,snr_idx, mod_idx) = mean(ber_1);
%         BER(2,snr_idx, mod_idx) = mean(ber_2);
        BER(3,snr_idx, mod_idx) = mean(ber_3);
        BER(4,snr_idx, mod_idx) = mean(ber_4);
        BER(5,snr_idx, mod_idx) = mean(ber_5);
        BER(6,snr_idx, mod_idx) = mean(ber_6);
        BER(7,snr_idx, mod_idx) = mean(ber_13);
        BER(8,snr_idx, mod_idx) = mean(ber_14);
        BER(9,snr_idx, mod_idx) = mean(ber_15);
        BER(10,snr_idx, mod_idx) = mean(ber_16);
    end
end




% Define colors for each modulation scheme
colours = [0.0000 0.4470 0.7410;
           0.8500 0.3250 0.0980;
           0.9290 0.6940 0.1250;
           0.4940 0.1840 0.5560;
           0.4660 0.6740 0.1880;
           0.3010 0.7450 0.9330;
           0.6350 0.0780 0.1840;
           0.0000 0.0000 0.0000];
% Create a new figure
figure;
savefig = 0;
linewidth = 2;
fontname = 'Times New Roman';
fontsize = 20;
legendsize = 14;
markersize = 10;
set(gcf, 'position', [0 0 800 600]);
root_save = 'C:\Program Files\UFRJ\TCC\images\new_model\';

% Plot BER vs SNR for each modulation scheme
for mod_idx = 1:num_mod
    % Plot BER for User 1 (Partial CSI-Test)
    semilogy(SNR_dB, BER(3, :, mod_idx), '-', 'color', colours(1, :), 'linewidth', linewidth, 'DisplayName', sprintf('User 1 Partial CSI'));
    hold on;
    % Plot BER for User 1 (Full CSI-Test)
    semilogy(SNR_dB, BER(4, :, mod_idx), '--', 'color', colours(1, :), 'linewidth', linewidth, 'DisplayName', sprintf('User 1 Full CSI'));
    % Plot BER for User 2 (Partial CSI-Test)
    semilogy(SNR_dB, BER(7, :, mod_idx), '-o', 'color', colours(2, :), 'linewidth', linewidth, 'MarkerSize', markersize, 'DisplayName', sprintf('User 2 Partial CSI'));
    % Plot BER for User 2 (Full CSI-Test)
    semilogy(SNR_dB, BER(8, :, mod_idx), '--o', 'color', colours(2, :), 'linewidth', linewidth, 'MarkerSize', markersize, 'DisplayName', sprintf('User 2 Full CSI'));
end

xlabel('SNR (in dB)', 'fontname', fontname, 'fontsize', fontsize);
ylabel('BER', 'fontname', fontname, 'fontsize', fontsize);
xlim([-10 35]);

legend('show', 'fontname', fontname, 'fontsize', legendsize, 'Location', 'southwest');
legend box off

set(gca, 'fontname', fontname, 'fontsize', fontsize);

% Set y-axis limits and ticks
%ylim([1e-3, 1]); % Set y-axis limits if needed
%yticks([1e-5, 1e-3, 1e-2, 1e-1, 1]); % Define specific tick marks

if savefig == 1
    saveas(gcf, [root_save 'new_model_BER_vs_SNR_dl_4QAM_40SNR_100epochs_test'], 'fig');
    saveas(gcf, [root_save 'new_model_BER_vs_SNR_dl_4QAM_40SNR_100epochs_test'], 'png');
    saveas(gcf, [root_save 'new_model_BER_vs_SNR_dl_4QAM_40SNR_100epochs_test'], 'epsc2');
end

save('ber_epochs_20', 'BER')
