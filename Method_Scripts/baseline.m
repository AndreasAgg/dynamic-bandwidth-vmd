function [u, u_hat, omega] = baseline(signal, K)
% Baseline Approach
% Code developer: Andreas Angelou
% 
% Detect frequeny peaks and apply bandpass filter. 
% 
% Implementation 
% --------------------------------
% We initialize the center frequencies to the frequency peaks of the signal 
% and we run the DB-VMD for only one iteration. 

%% Check inputs
arguments
    signal (1,:) double
    K double
end

%% Define parameters
tau_ab = 0.1;
tau_l = 0.1;
DC = 0;
init = 3;
tol = 1e-7;
max_it = 2;

%% Run DB-VMD for one iteration
[u, u_hat, omega] = DB_VMD_one(signal, tau_ab, tau_l, K, DC, init, tol, max_it);

end

%% DB-VMD for one iteration
function [u, u_hat, omega] = DB_VMD_one(signal, tau_ab, tau_l, K, DC, init, tol, max_it, opts)
% DB_VMD_one: DB-VMD one iteration (max_it = 2)
%% Check input
arguments
    signal (1,:) double
    tau_ab double
    tau_l double
    K double
    DC int8 
    init int8
    tol double
    max_it double
    opts.att double = -20
end

if ~all(opts.att < 0 & isa(opts.att, 'double'))
    error("att must be negative scaler")
end

%% Preparations
att_linear = db2mag(opts.att);

% Period and sampling frequency of input signal
N = length(signal);
Fs = N;

f_mirror = nan(2*N, 1);
f_mirror(1:N/2) = signal(N/2:-1:1);
f_mirror(N/2+1:3*N/2) = signal;
f_mirror(3*N/2+1:2*N) = signal(N:-1:N/2+1);
ext_signal = f_mirror;

% Time Domain 0 to T (of mirrored signal)
N_ext = length(ext_signal);
t = (1:N_ext)/N_ext;

% Spectral Domain discretization
freqs = t-0.5-1/N_ext; % [-0.5, ... , 0.495]

% Start with empty dual variables for bandwidth update
rho_alpha = zeros(max_it, 1);
rho_beta = zeros(max_it, 1);
rho_alpha(1) = 1;
rho_beta(1) = 16*K^2*(1/att_linear-1);

% Construct and center f_hat
f_hat = fftshift(fft(ext_signal));
f_hat_plus = f_hat;
f_hat_plus(1:N_ext/2) = 0; % Make negative freqs equal to zero
f_hat_plus = transpose(f_hat_plus);

% Matrix keeping track of every iterant
u_hat_plus = zeros(max_it, length(freqs), K);

% Initialization of omega_k at frequency peaks
omega_plus = nan(max_it, K);
omega_plus(1,:) = initializeCentralFreqs(init, K, f_hat_plus, freqs);

% If DC mode imposed, set its omega to 0
if DC
    omega_plus(1,1) = 0;
end

% Start with empty dual variables
lambda_hat = zeros(max_it, length(freqs));

% Other inits
uDiff = tol+eps; % update step
it = 1; % loop counter
sum_uk = 0; % accumulator

uDiffs = nan(max_it, 1);

% Occupied BW
ob = obw(signal, Fs);

% Values of alpha and beta
alpha = 0;
beta = ob/K;

%% Main loop for iterative updates
while ( uDiff > tol && it < max_it ) % not converged and below iterations limit

    % Update first mode accumulator
    k = 1;
    sum_uk = u_hat_plus(it,:,K) + sum_uk - u_hat_plus(it,:,1);
    
    % Update spectrum of first mode through Wiener filter of residuals
    u_hat_plus(it+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(it,:)/2) ./ (1 + (1 + rho_beta(it)-rho_alpha(it)) * (freqs - omega_plus(it,k)).^2);

    % Update first omega if not held at 0
    if ~DC
        omega_plus(it+1,k) = (freqs(N_ext/2+1:N_ext)*(abs(u_hat_plus(it+1, N_ext/2+1:N_ext, k)).^2)')/sum(abs(u_hat_plus(it+1,N_ext/2+1:N_ext,k)).^2);
    end
    
    % Update of any other mode
    for k=2:K
        
        % Accumulator
        sum_uk = u_hat_plus(it+1,:,k-1) + sum_uk - u_hat_plus(it,:,k);
        
        % Mode spectrum
        u_hat_plus(it+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(it,:)/2) ./ (1 + (1 + rho_beta(it)-rho_alpha(it)) * (freqs - omega_plus(it,k)).^2);
        
        % Center frequencies
        omega_plus(it+1,k) = (freqs(N_ext/2+1:N_ext)*(abs(u_hat_plus(it+1, N_ext/2+1:N_ext, k)).^2)')/sum(abs(u_hat_plus(it+1,N_ext/2+1:N_ext,k)).^2);
        
    end
    collective_BW = K * 2 * sqrt((1/att_linear - 1) / (1 + rho_beta(it)-rho_alpha(it))) * Fs;

    % Dual ascent
    lambda_hat(it+1,:) = lambda_hat(it,:) + tau_l*(sum(u_hat_plus(it+1,:,:),3) - f_hat_plus);
    
    % Bandwidth update dual variables
    rho_alpha(it+1) = rho_alpha(it) + tau_ab*(alpha - collective_BW);
    rho_beta(it+1) = rho_beta(it) + tau_ab*(collective_BW - beta);
    
    uDiffs(it) = uDiff;
    
    % Loop counter
    it = it+1;
    
    % Converged yet?
    uDiff = eps;
    for k=1:K
        uDiff = uDiff + 1/N_ext*(u_hat_plus(it,:,k)-u_hat_plus(it-1,:,k))*conj((u_hat_plus(it,:,k)-u_hat_plus(it-1,:,k)))';
    end
    uDiff = abs(uDiff);
    
end

%% Postprocessing and cleanup

% Discard empty space if converged early
omega = omega_plus(1:max_it,:);

% Signal reconstruction
u_hat = zeros(N_ext, K);
u_hat((N_ext/2+1):N_ext,:) = squeeze(u_hat_plus(max_it,(N_ext/2+1):N_ext,:));
u_hat((N_ext/2+1):-1:2,:) = squeeze(conj(u_hat_plus(max_it,(N_ext/2+1):N_ext,:)));
u_hat(1,:) = conj(u_hat(end,:));

u = zeros(K,length(t));

for k = 1:K
    u(k,:)=real(ifft(ifftshift(u_hat(:,k))));
end

% Remove mirror part
u = u(:,N_ext/4+1:3*N_ext/4);

% Recompute spectrum
clear u_hat;
for k = 1:K
    u_hat(:,k)=fftshift(fft(u(k,:)))';
end

end

%% Functions 
function centralFreqs = initializeCentralFreqs(init, K, x, f)
% Initialization of central frequencies
switch init
    case 3
        N_ext = length(f);
        x = abs(x(N_ext/2+1:end));
        f = f(N_ext/2+1:end);
        centralFreqs = initialCentralFreqByFindPeaks(x, f, K);
    otherwise
        error("Baseline approach: init must be 3.")
end
end

function centralFreq = initialCentralFreqByFindPeaks(x,f,K)
% Initialize central frequencies by finding the locations of signal peaks
% in frequency domain by using findpeaks function. The number of peaks is
% determined by NumIMFs.
BW = 2/(length(f)); % bandwidth of signal
minBWGapIndex = 2*BW/f(2);

x(x<mean(x)) = mean(x);
TF = islocalmax(x,'MinSeparation',minBWGapIndex);
pkst = x(TF);
locst = f(TF);
numpPeaks = length(pkst);

% Check for DC component
if x(1) >= x(2)
    pks = zeros(numpPeaks+1,1);
    locs = pks;
    pks(2:length(pkst)+1) = pkst;
    locs(2:length(pkst)+1) = locst;
    pks(1) = x(1);
    locs(1) = f(1);
else
    pks = zeros(numpPeaks,1);
    locs = pks;
    pks(1:length(pkst)) = pkst;
    locs(1:length(pkst)) = locst;
end   

[~,index] = sort(pks,'descend');
centralFreq = 0.5*rand(K,1);

% Check if the number of peaks is less than number of IMFs
if length(locs) < K
    centralFreq(1:length(locs(index))) = locs;
else
    centralFreq(1:K) = locs(index(1:K));
end

centralFreq = sort(centralFreq);

end
