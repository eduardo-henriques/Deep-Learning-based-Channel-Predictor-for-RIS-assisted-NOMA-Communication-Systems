data = load('dataset_1.mat');
Y = data.Y;
K = size(Y,3);
    % Assuming Y is of size [2*N_t*N_r, M, K, SNR] 
    % and the first half of the first dimension contains the amplitudes
    
    SNR_values = size(Y, 4); % Get the number of SNR values
    time_steps = 1:K; % Prepare the time step axis
    
    figure; % Create a new figure for the plots
    
%    for snr = 1:SNR_values
        % Extract the amplitude of the received signals for both users
        amplitude_user1 = squeeze(Y(1, 1, :, 9, 1)).^2; % Amplitudes for user 1
        amplitude_user2 = squeeze(Y(1, 2, :, 9, 1)).^2; % Amplitudes for user 2

        % Normalize the amplitudes
        max_amplitude = max(amplitude_user1(:));
        normalized_amplitude_user1 = amplitude_user1 / max_amplitude;
        normalized_amplitude_user2 = amplitude_user2 / max_amplitude;
        
        savefig = 0;
        root_save = 'C:\Program Files\UFRJ\TCC\images\';
        linewidth  = 2;
        markersize = 10;
        fontname   = 'Times New Roman';
        fontsize   = 20;
        
        legend_alg = {'User 1', 'User 2'};
        
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
        
        set(gcf,'position',[0 0 800 600]);
        semilogy(time_steps, normalized_amplitude_user1, 'r-', 'color',colours(1,:),'linewidth',linewidth);
        hold on;
        semilogy(time_steps, normalized_amplitude_user2, 'b-','color',colours(2,:),'linewidth',linewidth);
        
        ylabel('Normalized Received Power','fontname',fontname,'fontsize',fontsize);
        xlabel('Time (in samples)','fontname',fontname,'fontsize',fontsize);
        
        legend(legend_alg,'fontname',fontname,'fontsize',fontsize,'location',location_4);
        legend boxoff;
        set(gca,'fontname',fontname,'fontsize',fontsize);

        if savefig == 1
            saveas(gcf,[root_save 'normalized_received_power_vs_timesteps_0db'], 'fig');
            saveas(gcf,[root_save 'normalized_received_power_vs_timesteps_0db'], 'png');
            saveas(gcf,[root_save 'normalized_received_power_vs_timesteps_0db'], 'epsc2');
        end
%    end
