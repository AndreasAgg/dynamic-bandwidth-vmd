%% Add "Method_Scripts" path
% Method_Scripts is the folder where DB-VMD and VMD are implemeneted 

init_pwd = pwd;
cd ..
addpath('Method_Scripts') 
cd(init_pwd)


%% Preparation 
clear;
close all;
clc;

% Time Domain 0 to N
N = 800;
Fs = N;
t = (1:N)/N;

start_freq = 1;
step = 3;
f1_arr = start_freq : step : N/2;

SR_arr_DB_VMD = ones(length(f1_arr), length(f1_arr), 3);
SR_arr_VMD = ones(length(f1_arr), length(f1_arr), 3);

rho_arr = [1/4, 1, 4];

alpha = 1000;
K = 2;
tau_ab = 0.1;
tau_l = 0.1;
DC = 0;
init = 3;
tol = 1e-7;
%% Calculation of error array DB-VMD
f1_it = 1;
for v_1 = f1_arr
    f2_it = 1;
    for v_2 = start_freq : step : v_1
        rho_it = 1;
        for rho = rho_arr
            a1 = randi(3)*rand(1,1);
            a2 = rho * a1;
            x_1 = a1 * (cos(2*pi*v_1*t));
            x_2 = a2 * (cos(2*pi*v_2*t));
            fsub = {};
            fsub{1} = x_1;
            fsub{2} = x_2;
            x = x_1 + x_2;

            [u, ~, omega] = DB_VMD(x, tau_ab, tau_l, K, DC, init, tol);
            [~, sortIndex] = sort(omega(end,:), 'descend');
            u = u(sortIndex,:);
            corr_arr = nan(K,1);
            for k=1:K
                corr_arr(k) = abs(xcorr(fsub{k},u(k,:),0,'normalized'));
            end
            SR_arr_DB_VMD(end-f2_it+1, f1_it, rho_it) = mean(corr_arr);
            
            [u, ~, omega] = VMD(x, alpha, tau_l, K, DC, init, tol);
            [~, sortIndex] = sort(omega(end,:), 'descend');
            u = u(sortIndex,:);
            corr_arr = nan(K,1);
            for k=1:K
                corr_arr(k) = abs(xcorr(fsub{k},u(k,:),0,'normalized'));
            end
            SR_arr_VMD(end-f2_it+1, f1_it, rho_it) = mean(corr_arr);

            rho_it = rho_it + 1;
        end
        f2_it = f2_it  + 1;
    end
    f1_it = f1_it + 1;
end
%% Visualization

for i=1:length(rho_arr)
    figure("Name", 'DB-VMD: rho=' + sprintf("%s", num2str(rho_arr(i)))); 
    imshow(SR_arr_DB_VMD(:,:,i), 'InitialMagnification', 'fit')
    title("DB-VMD", "Interpreter", "latex", 'FontSize', 30)
    xlabel("$v_1$", "Interpreter", "latex", 'FontSize', 20)
    ylabel("$v_2 < v_1$", "Interpreter", "latex", 'FontSize', 20)
    colormap gray
    colorbar
    caxis([0.3, 1]);
    colorbar off
    
    figure("Name", 'VMD: rho=' + sprintf("%s", num2str(rho_arr(i)))); 
    imshow(SR_arr_VMD(:,:,i), 'InitialMagnification', 'fit')
    title("VMD", "Interpreter", "latex", 'FontSize', 30)
    xlabel("$v_1$", "Interpreter", "latex", 'FontSize', 20)
    ylabel("$v_2 < v_1$", "Interpreter", "latex", 'FontSize', 20)
    colormap gray
    caxis([0.3, 1]);
    colorbar off

end