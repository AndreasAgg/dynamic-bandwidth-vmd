%% Add "Method_Scripts" path
% Method_Scripts is the folder where DB-VMD and VMD are implemeneted 

init_pwd = pwd;
cd ..
addpath('Method_Scripts') 
cd(init_pwd)

%% Preparations and parameters definition

clear; clc; close all

snr_arr = 20:-2:10;
tau_l_arr = [0.1, 0.5];

max_it = 500; % Number of iterations

% Parameters
alpha = 1000;   % VMD bandwidth factor
K = 3;          % Components' count
tau_ab = 0.1;   % Bandwidth rate of change (DB-VMD)
DC = 0;         % DC impose (0 for none)
init = 3;       % Central frequencies initialization
tol = 1e-7;     % Stopping criteria tolerance


N = 500; % Signal length
n = (1:N)';

% Generating Hanning windows
L_arr = [500, 125, 100];
d_arr = [250, 125, 375];
hann_windows = nan(N,K);

for i=1:K
    temp = zeros(N,1);
    low = d_arr(i) - L_arr(i)/2 + 1;
    high = d_arr(i) + L_arr(i)/2;
    temp(low: high) = hann(L_arr(i));
    hann_windows(:, i) = temp;
end

% Success rate vector for max_it iterations
SR_DB_vmd_vec = nan(max_it, 1);
SR_vmd_vec = nan(max_it, 1);

% Success rate array for (SNR, tau_l) pairs 
SR_DB_vmd_arr = nan(length(snr_arr), length(tau_l_arr));
SR_vmd_arr = nan(length(snr_arr), length(tau_l_arr));

%% Noise robustness experiment
for i_tau_l = 1:length(tau_l_arr)
    tau_l = tau_l_arr(i_tau_l);
    
    for i_snr = 1:length(snr_arr)
        for it = 1:max_it
            
            % Printing progress
            if mod(it, 50) == 0
                fprintf("iteration: %d/%d - SNR: %.2f (%d/%d) - tau_l: %.2f (%d/%d)\n", ...
                    it, max_it, ...
                    snr_arr(i_snr), i_snr, length(snr_arr), ...
                    tau_l_arr(i_tau_l), i_tau_l, length(tau_l_arr))
            end            
            
            % Signal generation
            omega_arr = unifrnd(0,pi,3,1);
            A_arr = unifrnd(0.5,1.5,3,1);
            fsub = cell(K,1);
            s = zeros(N,1);
            for i=1:K
                fsub{i} = hann_windows(:, i) .* A_arr(i) .* cos(omega_arr(i).*n);
                s = s + hann_windows(:, i) .* A_arr(i) .* cos(omega_arr(i).*n);
            end
            [~, sortIndex] = sort(omega_arr);
            fsub = fsub(sortIndex);
            s = awgn(s, snr_arr(i_snr),"measured"); % Add noise
            
            % DB-VMD applied
            [u, ~, omega] = DB_VMD(s, tau_ab, tau_l, K, DC, init, tol);
            [~, sortIndex] = sort(omega(end,:));
            u = u(sortIndex,:);

            % DB-VMD Success rate 
            corr_arr = nan(K,1);
            for k=1:K
                corr_arr(k) = xcorr(fsub{k},u(k,:),0,'normalized');
            end
            SR_DB_vmd_vec(it) = mean(corr_arr);


            % VMD applied
            [u, ~, omega] = VMD(s, alpha, tau_l, K, DC, init, tol);
            [~, sortIndex] = sort(omega(end,:));
            u = u(sortIndex,:);

            % VMD Success rate 
            corr_arr = nan(K,1);
            for k=1:K
                corr_arr(k) = xcorr(fsub{k},u(k,:),0,'normalized');
            end
            SR_vmd_vec(it) = mean(corr_arr);
        end
        
        % Average SR over max_it success rates
        SR_DB_vmd_arr(i_snr, i_tau_l) = mean(SR_DB_vmd_vec);
        SR_vmd_arr(i_snr, i_tau_l) = mean(SR_vmd_vec);
    end  
end
%% Results 

for i=1:length(tau_l_arr)
    figure("Name", sprintf("tau_l = %.2f", tau_l_arr(i)))
    semilogx(snr_arr, SR_DB_vmd_arr(:, i), 'k-')
    hold on
    semilogx(snr_arr, SR_vmd_arr(:, i), 'k:')
    xlabel("SNR", 'FontSize', 15)
    ylabel("Success Rate", 'FontSize', 20)
    legend("DB-VMD","VMD", "Location", "Best", 'FontSize',12)
    set(gca, 'XGrid', 'on', 'XMinorGrid', 'on');
    title("Success Rate vs SNR", 'FontSize', 20)
end


