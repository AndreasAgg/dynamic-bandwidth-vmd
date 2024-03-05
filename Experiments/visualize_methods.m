%% Script for visualizing DB-VMD and/or VMD
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

% Set the method variable below to "DB-VMD" or "VMD" 
method = "DB-VMD";
if ~all(strcmp(method, "VMD") | strcmp(method, "DB-VMD"))
    error("method must be VMD or DB-VMD")
end

% Time Domain 0 to 1
N = 1000;
Ts = 1/N;
t = (1:N)/N;
Fs = N;
freqs = (t-0.5-Ts)/Ts;

% Center frequencies of components
f_1 = 100;
f_2 = 200;
f_3 = 400;

% Components
v_1 = cos(2*pi*f_1*t);
v_2 = 1/4*(cos(2*pi*f_2*t));
v_3 = 1/16*(cos(2*pi*f_3*t));

% For visualization purposes
fsub = {};
wsub = {};
fsub{1} = v_1;
fsub{2} = v_2;
fsub{3} = v_3;

wsub{1} = 2*pi*f_1;
wsub{2} = 2*pi*f_2;
wsub{3} = 2*pi*f_3;

% Pure signal and noise
x = v_1 + v_2 + v_3;
noise = 0.05*randn(size(x));

% Composite signal, including noise
x = x + noise; 

% FFT of composite signal
F = fftshift((fft(x)));

% Parameters for both methods
K = 3;
tau_l = 0.1;
DC = 0;
init = 3;
tol = 1e-6;

if strcmp(method, "DB-VMD")
    % BW rate of change for DB-VMD
    tau_ab = 0.1;
else
    % BW factor for VMD
    alpha = 1000;
end
%% Run DB-VMD (or VMD) code and visualize progress
if strcmp(method, "DB-VMD")
    [u, u_hat, omega] = DB_VMD(x, tau_ab, tau_l, K, DC, init, tol, "viz_progress", 1);
else
    % VMD method
    [u, u_hat, omega] = VMD(x, alpha, tau_l, K, DC, init, tol, "viz_progress", 1);
end
%% Visualization
% For convenience here: Order omegas increasingly and reindex u and u_hat
[~, sortIndex] = sort(omega(end,:));
omega = omega(:,sortIndex);
u_hat = u_hat(:,sortIndex);
u = u(sortIndex,:);
linestyles = {'b', 'g', 'm', 'c', 'y', 'r', 'k'};

% Plot composite signal
figure('Name', method + ": " + 'Composite input signal' );
plot(t, x, 'k');
title(method + ": " + "Composite Signal")
set(gca, 'XLim', [0 1]);
xlabel("Time")

% Plot FFT of composite signal
figure('Name', method + ": " + 'Input signal spectrum' );
loglog(2*pi*freqs(N/2+1:end), abs(F(N/2+1:end)), 'k');
set(gca, 'XLim', [1 N/2*2*pi], 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylims = get(gca, 'YLim');
hold on;
for sub = 1:length(wsub)
    xline(wsub{sub}, 'k--')
end
set(gca, 'YLim', ylims);
xlabel("Cyclic frequencies $\omega$", "interpreter", "latex")
title(method + ": " + "Input signal spectrum")

% Plot evolution of center frequencies
figure('Name', method + ": " + 'Evolution of center frequencies omega');
for k=1:K
    semilogx(2*pi/Ts*omega(:,k), 1:size(omega,1), 'Color', linestyles{k});
    hold on;
end
set(gca, 'YLim', [1,size(omega,1)]);
set(gca, 'XLim', [2*pi,0.5*2*pi/Ts], 'XGrid', 'on', 'XMinorGrid', 'on');
title(method + ": " + "Evolution of center frequencies omega")
xlabel("Central frequencies $\omega_k$", "interpreter", "latex")
ylabel("Iteration count")

% Plot input spectrum and modes' spectra
figure('Name', method + ": " + 'Spectral decomposition');
loglog(2*pi*freqs(N/2+1:end), abs(F(N/2+1:end)), 'k:');
set(gca, 'XLim', [1 N/2]*pi*2, 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off');
hold on;
for k = 1:K
    loglog(2*pi*freqs(N/2+1:end), abs(u_hat(N/2+1:end,k)), 'Color', linestyles{k});
end
set(gca, 'YLim', ylims);
title(method + ": " + "Spectral decomposition")

% Plot reconstructed modes and true components
for k = 1:K
    figure('Name', method + ": " + ['Reconstructed mode ' num2str(k)]);
    plot(t,u(k,:), 'Color', linestyles{k});   hold on;
    if ~isempty(fsub)
        plot(t, fsub{min(k,length(fsub))}, 'k:');
    end
    set(gca, 'XLim', [0 1]);
    title(method + ": " + ['Reconstructed mode ' num2str(k)])
    
end

% Plot pure signal and the supercomposed signal
figure('Name', method + ": " + 'Comparison: Modes - Real Signal')
plot(t, x - noise, 'k:', 'Linewidth', 1.5)
hold on
plot(t, sum(u), 'r')
set(gca, 'XLim', [0 1], 'XGrid', 'on', 'YGrid', 'on');
xlabel("Time (t)")

% Calculate the mean normalized cross-correlation
corr_arr = nan(K,1);
for k=1:K
    corr_arr(k) = xcorr(fsub{k},u(k,:),0,'normalized');
end
SR = mean(corr_arr);
title(method + ": " + sprintf("Mean cross-correlation: %s", string(SR)))
legend("Noise-free input signal", "Reconstructed signal")
