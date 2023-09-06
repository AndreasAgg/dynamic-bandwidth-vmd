function [u, u_hat, omega] = VMD(signal, alpha, tau, K, DC, init, tol, opts)
% Variational Mode Decomposition
% Authors: Konstantin Dragomiretskiy and Dominique Zosso
% zosso@math.ucla.edu --- http://www.math.ucla.edu/~zosso
% Initial release 2013-12-12 (c) 2013
%
% Input and Parameters:
% ---------------------
% signal  - the time domain signal (1D) to be decomposed
% alpha   - the balancing parameter of the data-fidelity constraint
% tau     - time-step of the dual ascent ( pick 0 for noise-slack )
% K       - the number of modes to be recovered
% DC      - true if the first mode is put and kept at DC (0-freq)
% init    - 0 = all omegas start at 0
%           1 = all omegas start uniformly distributed
%           2 = all omegas initialized randomly
%           3 = all omegas initialized as the peak locations of the signal in the frequency domain
% tol     - tolerance of convergence criterion; typically around 1e-6
% opts    - Options used with the "'Name', Value" format
%               VMD(..., 'viz_end', Value1, 'viz_progress', Value2)
%                   viz_progress: visualization while the algorithm is in progress
%                   viz_end: visualization after the end of the algorithm 
%
% Output:
% -------
% u       - the collection of decomposed modes
% u_hat   - spectra of the modes
% omega   - estimated mode center-frequencies
%
% When using this code, please do cite our paper:
% -----------------------------------------------
% K. Dragomiretskiy, D. Zosso, Variational Mode Decomposition, IEEE Trans.
% on Signal Processing (in press)
% please check here for update reference: 
%          http://dx.doi.org/10.1109/TSP.2013.2288675
%% Check inputs
arguments
    signal (1,:) double
    alpha double
    tau double
    K int8
    DC int8 
    init int8
    tol double
    opts.viz_end int8 = 0
    opts.viz_progress int8 = 0
end


if ~all(opts.viz_end == 0 | opts.viz_end == 1)
    error("viz_end must be 0 or 1")
end

if ~all(opts.viz_progress == 0 | opts.viz_progress == 1)
    error("viz_progress must be 0 or 1")
end

%% ---------- Preparations
% Period and sampling frequency of input signal
N = length(signal);

% Extend the signal by mirroring
f_mirror(1:N/2) = signal(N/2:-1:1);
f_mirror(N/2+1:3*N/2) = signal;
f_mirror(3*N/2+1:2*N) = signal(N:-1:N/2+1);
ext_signal = f_mirror;

% Time Domain 0 to T (of mirrored signal)
N_ext = length(ext_signal);
t = (1:N_ext)/N_ext;

% Spectral Domain discretization
freqs = t-0.5-1/N_ext; % [-0.5, ... , 0.495]

% Maximum number of iterations (if not converged yet, then it won't anyway)
max_it = 1000;

% For future generalizations: individual alpha for each mode
Alpha = alpha*ones(1,K);

% Construct and center f_hat
f_hat = fftshift((fft(ext_signal)));
f_hat_plus = f_hat;
f_hat_plus(1:N_ext/2) = 0; % Make negative freqs equal to zero

% Matrix keeping track of every iterant // could be discarded for mem
u_hat_plus = zeros(max_it, length(freqs), K);

% Initialization of omega_k
omega_plus = nan(max_it, K);
omega_plus(1,:) = initializeCentralFreqs(init, K, f_hat_plus, freqs);

% if DC mode imposed, set its omega to 0
if DC
    omega_plus(1,1) = 0;
end

% start with empty dual variables
lambda_hat = zeros(max_it, length(freqs));

% other inits
uDiff = tol+eps; % update step
it = 1; % loop counter
sum_uk = 0; % accumulator

%% Viz figure
method = "VMD";
if opts.viz_progress || opts.viz_end
    figure('Name', 'Visualization ' + method)
end

%% ----------- Main loop for iterative updates
while ( uDiff > tol &&  it < max_it ) % not converged and below iterations limit
    
    % Update first mode accumulator
    k = 1;
    sum_uk = u_hat_plus(it,:,K) + sum_uk - u_hat_plus(it,:,1);
    
    % Update spectrum of first mode through Wiener filter of residuals
    u_hat_plus(it+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(it,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(it,k)).^2);
    
    % Update first omega if not held at 0
    if ~DC
        omega_plus(it+1,k) = (freqs(N_ext/2+1:N_ext)*(abs(u_hat_plus(it+1, N_ext/2+1:N_ext, k)).^2)')/sum(abs(u_hat_plus(it+1,N_ext/2+1:N_ext,k)).^2);
    end
    
    % Update of any other mode
    for k=2:K
        
        % Accumulator
        sum_uk = u_hat_plus(it+1,:,k-1) + sum_uk - u_hat_plus(it,:,k);
        
        % Mode spectrum
        u_hat_plus(it+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(it,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(it,k)).^2);
        
        % Center frequencies
        omega_plus(it+1,k) = (freqs(N_ext/2+1:N_ext)*(abs(u_hat_plus(it+1, N_ext/2+1:N_ext, k)).^2)')/sum(abs(u_hat_plus(it+1,N_ext/2+1:N_ext,k)).^2);
        
    end
    
    % Dual ascent
    lambda_hat(it+1,:) = lambda_hat(it,:) + tau*(sum(u_hat_plus(it+1,:,:),3) - f_hat_plus);
    
    if opts.viz_progress
        visualizationFunction(omega_plus(it+1, :), freqs, f_hat_plus, alpha, "VMD");
    end
    
    % Loop counter
    it = it+1;
        
    % Converged yet?
    uDiff = eps;
    for i=1:K
        uDiff = uDiff + 1/N_ext*(u_hat_plus(it,:,i)-u_hat_plus(it-1,:,i))*conj((u_hat_plus(it,:,i)-u_hat_plus(it-1,:,i)))';
    end
    uDiff = abs(uDiff);
    
end
it = it - 1;

% Visualize end
if opts.viz_end
    visualizationFunction(omega_plus(it, :), freqs, f_hat_plus, alpha, "VMD");
    fprintf("VMD iterations: %d\n", it)
end

%% ------ Postprocessing and cleanup

% Discard empty space if converged early
max_it = min(max_it,it);
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
    case 0
        centralFreqs = initialCentralFreqToZero(K);
    case 1
        centralFreqs = initialCentralFreqUniformally(K);
    case 2
        Ts = 1 / length(x);
        centralFreqs = initialCentralFreqRandomly(K, Ts);
    case 3
        N_ext = length(f);
        x = abs(x(N_ext/2+1:end));
        f = f(N_ext/2+1:end);
        centralFreqs = initialCentralFreqByFindPeaks(x, f, K);
    otherwise
        error("init is not between 0 and 3")
end
end


function centralFreqs = initialCentralFreqToZero(K)

centralFreqs = zeros(1, K);

end

function centralFreqs = initialCentralFreqUniformally(K)

centralFreqs = nan(K, 1);

for i = 1:K
    centralFreqs(i) = (0.5/K)*(i-1);
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

function centralFreqs = initialCentralFreqRandomly(K, Ts)

centralFreqs = sort(exp(log(Ts) + (log(0.5)-log(Ts))*rand(1,K)));

end
