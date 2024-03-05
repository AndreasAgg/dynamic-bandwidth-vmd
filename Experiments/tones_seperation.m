%% Tones separation script for Dynamic Bandwidth VMD
% -------------------------------------------------------------------------
% Written by:
% Andreas G. Angelou, aangelou@auth.gr
% Tested with MATLAB R2020b
% -------------------------------------------------------------------------
%% Add "Method_Scripts" path
% Method_Scripts is the folder where DB-VMD and VMD are implemeneted 

init_pwd = pwd;
cd ..
addpath('Method_Scripts') 
cd(init_pwd)

%% Preparation and parameters definition
clear; close all; clc;

N = 800; % Signal length
t = (1:N)/N; % Time vector

% Frequencies array
start_freq = 1;
step = 3;
f1_arr = start_freq : step : N/2;

% Initialization of Success Rate arrays
SR_arr_DB_VMD = ones(length(f1_arr), length(f1_arr), 3);
SR_arr_VMD = ones(length(f1_arr), length(f1_arr), 3);

% Rho array
rho_arr = [1/4, 1, 4];

% Parameters
alpha = 1000;   % VMD bandwidth factor
K = 2;          % Components' count
tau_ab = 0.1;   % Bandwidth rate of change (DB-VMD)
tau_l = 0.1;    % Data-fidelity factor
DC = 0;         % DC impose (0 for none)
init = 3;       % Central frequencies initialization
tol = 1e-7;     % Stopping criteria tolerance
%% Tone seperation experiment
f1_it = 1;
for v_1 = f1_arr
    
    % Printing progress
    fprintf("Progress: %d/%d\n", f1_it, length(f1_arr));
    
    f2_it = 1;
    for v_2 = start_freq : step : v_1
        rho_it = 1;
        for rho = rho_arr
            
            % Signal generation
            a1 = randi(3) * rand(1, 1);
            a2 = rho * a1;
            x_1 = a1 * (cos(2 * pi * v_1 * t));
            x_2 = a2 * (cos(2 * pi * v_2 * t));
            fsub = {};
            fsub{1} = x_1;
            fsub{2} = x_2;
            x = x_1 + x_2;

            % DB-VMD applied
            [u, ~, omega] = DB_VMD(x, tau_ab, tau_l, K, DC, init, tol);
            [~, sortIndex] = sort(omega(end, :), 'descend');
            u = u(sortIndex, :);
            
            % DB-VMD Success rate 
            corr_arr = nan(K, 1);
            for k=1:K
                corr_arr(k) = abs(xcorr(fsub{k}, u(k,:), 0, 'normalized'));
            end
            SR_arr_DB_VMD(end - f2_it + 1, f1_it, rho_it) = mean(corr_arr);
            
            
            % VMD applied
            [u, ~, omega] = VMD(x, alpha, tau_l, K, DC, init, tol);
            [~, sortIndex] = sort(omega(end, :), 'descend');
            u = u(sortIndex, :);
            
            % VMD Success rate 
            corr_arr = nan(K, 1);
            for k=1:K
                corr_arr(k) = abs(xcorr(fsub{k}, u(k,:), 0, 'normalized'));
            end
            SR_arr_VMD(end - f2_it + 1, f1_it, rho_it) = mean(corr_arr);

            rho_it = rho_it + 1;
        end
        f2_it = f2_it  + 1;
    end
    f1_it = f1_it + 1;
end
%% Results

for i=1:length(rho_arr)
    figure("Name", 'DB-VMD: rho=' + sprintf("%s", num2str(rho_arr(i)))); 
    imshow(SR_arr_DB_VMD(:, :, i), 'InitialMagnification', 'fit')
    title("DB-VMD", "Interpreter", "latex", 'FontSize', 30)
    xlabel("$v_1$", "Interpreter", "latex", 'FontSize', 20)
    ylabel("$v_2 < v_1$", "Interpreter", "latex", 'FontSize', 20)
    colormap hot
    caxis([0.3, 1]);
    colorbar
    
    figure("Name", 'VMD: rho=' + sprintf("%s", num2str(rho_arr(i)))); 
    imshow(SR_arr_VMD(:, :, i), 'InitialMagnification', 'fit')
    title("VMD", "Interpreter", "latex", 'FontSize', 30)
    xlabel("$v_1$", "Interpreter", "latex", 'FontSize', 20)
    ylabel("$v_2 < v_1$", "Interpreter", "latex", 'FontSize', 20)
    colormap hot
    caxis([0.3, 1]);
    colorbar 

end