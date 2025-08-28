% Notes
% - The main goal is to implement a digital-to-analog converter
% (DAC) that's able to reconstruct the selected signals into one received
% signal in time.
% - Parameters were changed with respect to the last implementations,
% notably the correction of the carrier frequency (to match UWB standards)
% and the sampling frequency (to obey Nyquist's theorem).

% Functions ===============================================================

% Generate a single Gaussian Pulse (p(t))
function p = GaussianPulse(t, Vg, sig)
    % Signal
    p = Vg*exp(-(t.^2) / (2 * sig^2));
end

% Received RF pulse signal
function [t_rel, g_prime] = Received(t, Tp, M, fc, Vg, sig, alpha, c, R)

    % Parameters
    g_prime = zeros(size(t));
    Ts = t(2) - t(1); % assume uniform sampling
    N = round(Tp / Ts); % samples per pulse

    for m_idx = 1:M
        % Time window for single pulse
        t_start = (m_idx-1)*Tp;
        t_window = t - t_start;

        % Center of the fast-time window
        t_center = (N-1)/2 * Ts;

        % Target delay at pulse center
        t_abs_center = t_start + t_center;
        x_center = 2*sin(1*t_abs_center); % motion of the chest
        tau_center = 2*(R + x_center)/c;

        % Fast-time relative to pulse center
        t_rel = t_window - t_center - tau_center;

        % Gaussian pulse
        p_rf = GaussianPulse(t_rel, Vg, sig);

        % Carrier modulation at relative fast-time
        g_prime = g_prime + alpha * p_rf .* cos(2*pi*fc*t_rel);
    end
end

% Single-point model
function [gs_prime, r] = ADC(Tp, M, fc, Vg, sig, alpha, c, R, fs, N)

    % Matrix with dimensions pulses x samples
    gs_prime = zeros(M, N);

    for m_idx = 1:M
        % Pulse start time
        t_start = (m_idx-1)*Tp;

        % Target motion (example sinusoidal)
        x_center = 20*sin(1*t_start);
        tau_center = 2*(R + x_center)/c; % round-trip delay

        for n_idx = 1:N
            % Absolute fast-time
            t_fast = (n_idx-1)/fs;

            % Relative time wrt the delay
            t_rel = t_fast - tau_center;

            % Gaussian pulse
            p_mn = GaussianPulse(t_rel, Vg, sig);

            % Carrier modulation at relative fast-time
            gs_prime(m_idx, n_idx) = alpha * p_mn * cos(2*pi*fc*t_rel);
        end
    end

    % Marking the correct range bin
    r = 2*R*fs/c;
end

% DAC conversion
function [t_dac, g_dac, g_dac_matrix] = DAC(gs_prime, fs, up_factor)

    [M, N] = size(gs_prime);

    % Flatten pulses to chronological 1-D sequence
    g1d = gs_prime.'; % N x M
    g1d = g1d(:); % (M*N) x 1

    % Resample once using sinc-like interpolation
    g_up = resample(g1d, up_factor,1);  % column vector

    % Time vector for reconstructed waveform
    fs_up = fs * up_factor;
    t_dac = (0:length(g_up)-1).' / fs_up;

    % Visualization of per-sample pulses
    N_up = N * up_factor;
    if mod(length(g_up), N_up) == 0 && M*N_up == length(g_up)
        g_dac_matrix = reshape(g_up, N_up, M).';
    else
        g_dac_matrix = [];
    end

    g_dac = g_up;

end


% Variables ===============================================================

[N, M] = deal(400, 500); % Number of pulses, Total number of samples in each pulse interval
fs = 10e9; % Sampling frequency (in Hz)
Ts = 1/fs; % Sampling period (in sec)
fc = 5e9; % Carrier frequency (in Hz)
Tp = N/fs; % Pulse repetition interval (PRI, in sec)
t = 0:1/fs:M*Tp; % Time span, arbitrary
fast_time = (0:N-1)/fs * 1e6; % (in msec)
slow_time = (0:M-1) * Tp; % (in sec)
Vg = 1; % Amplitude
sig = 1e-9; % Standard deviation
c = physconst('LightSpeed'); % Speed of ligth (in m/s)
R = 2; % Distance to the target
alpha = 1; % Channel gain
fp = fs/100; % Passband frequency, arbitrary == only produces relevant results for some specific values
SNR_dB = 50; % Signal-to-noise ratio (in dB)
range_bins = [10, 50, 100, 150, 200, 250]; % Arbitrary
up_factor = 10; % DAC parameter

% Main ====================================================================

% g^p(t)
[t_rel, g_prime] = Received(t, Tp, M, fc, Vg, sig, alpha, c, R);

% gs_prime
[gs_prime, r] = ADC(Tp, M, fc, Vg, sig, alpha, c, R, fs, N);

% Reconstruct waveform
[t_dac, g_dac, g_dac_matrix] = DAC(gs_prime, fs, up_factor);

% Plots ===================================================================

% g_prime(t)
figure;
plot(t*1e3, g_prime);
xlabel('Time [ms]'); ylabel('Amplitude');
title('Received signal g^p(t)');

%gs_prime(m,n)
figure;
imagesc(1:N, 1:M, abs(gs_prime));
xline(r, 'LineWidth', 3);
xlabel('Fast time index (n)');
ylabel('Slow time index (m)');
title('ADC |g_{s}^{p}(m,n)|');
colorbar;

% DAC
figure;
plot(t_dac*1e6, g_dac);
xlabel('Time [\mus]'); ylabel('Amplitude');
title('Samples of pulses reconstructed by the DAC');

% Testing =================================================================

%Choosing a representative pulse index
m_idx = 10;
pulse = gs_prime(m_idx,:);

% Raw fast-time samples for the chosen pulse
figure; subplot(2,1,1);
plot((0:length(pulse)-1)/fs*1e6, pulse);
xlabel('fast-time (\mus)'); ylabel('amplitude'); title(['Raw pulse m=' num2str(m_idx)]);

% Analytic envelope to check if it's symmetric
env = abs(hilbert(pulse));
subplot(2,1,2);
plot((0:length(pulse)-1)/fs*1e6, env);
xlabel('fast-time (\mus)'); ylabel('envelope'); title('Envelope (abs(hilbert))');
