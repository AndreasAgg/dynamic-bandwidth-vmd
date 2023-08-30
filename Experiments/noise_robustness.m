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

max_it = 500;
K = 3;
tau_ab = 0.1;
DC = 0;
init = 3;
tol = 1e-7; 
alpha = 1000;


N = 500; % Signal length
n = (1:N)';
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

%% Actual Experiment
for i_tau_l = 1:length(tau_l_arr)
    tau_l = tau_l_arr(i_tau_l);
    
    for i_snr = 1:length(snr_arr)
        for it = 1:max_it
            omega_arr = unifrnd(0,pi,3,1);
            A_arr = unifrnd(0.5,1.5,3,1);

            % Signal generation
            fsub = cell(K,1);
            s = zeros(N,1);
            for i=1:K
                fsub{i} = hann_windows(:, i) .* A_arr(i) .* cos(omega_arr(i).*n);
                s = s + hann_windows(:, i) .* A_arr(i) .* cos(omega_arr(i).*n);
            end
            [~, sortIndex] = sort(omega_arr);
            fsub = fsub(sortIndex);

            s = awgn(s, snr_arr(i_snr),"measured");
            
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
            
            % Printing progress
            if mod(it, 50) == 0
                fprintf("iteration: %d - SNR: %.2f - tau_l: %.2f\n", it, snr_arr(i_snr), tau_l_arr(i_tau_l))
            end
        end
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


