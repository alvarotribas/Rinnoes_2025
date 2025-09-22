% Notes:
% - For simplicity, no range bin is being selected. The multi-point model
% contains all the bins from all the points it considers. This may lead to
% errors, which are guessed to be minor, compared to the strength of the
% actual received signal.

% Functions ===============================================================
% Generate a single Gaussian Pulse (p(t))
function p = GaussianPulse(t, Vg, sig)
    % Signal
    p = Vg*exp(-(t.^2) / (2 * sig^2));
end

% Complex carrier signal
function [n,s] = ComplexCarrier(Vg, N, fc, fs)
    n = 1:N;
    s = Vg * exp(-1i*2*pi*fc*(n/fs));
end

% Multi-point model
function [m, n, q, gq_prime, yH_mn, yH_mn_unf] = ADC_MultiTarget(M, N, Tp, fc, fp, fs, Vg, sig, c, SNR_dB)
    Ts = 1/fs;
    gq_prime = zeros(M, N);

    % Multi-target config
    Q = 50; % Number of targets
    chest = 5; % Position of the chest wrt the radar
    R_all = linspace(0.6, 4, Q); % Ranges of the human body, linear from closest to furthest
    alpha_all = normpdf(linspace(0.3, 0.9, Q), chest/10, 0.5); % Channel gains, linear within arbitrary range
    a_all = normpdf(linspace(1, 10, Q), chest, 2); % Movement constant, normal distributed around torax region
    b_all = normpdf(linspace(0.1, 1, Q), chest/10, 0.2); % Movement frequency, also normal around torax
    
    % Determining p_q(m,n)
    for q = 1:Q
        for m_idx = 1:M
            for n_idx = 1:N
                % Parameters of each received pulse
                t_mn = (m_idx-1)*Tp + (n_idx-1)*Ts;
                R_q = R_all(q);
                alpha_q = alpha_all(q);
                a = a_all(q);
                b = b_all(q);

                % Body movement and propagation delay
                x_q = a*sin(b * t_mn);
                tau_q = 2 * (R_q + x_q) / c;
                
                % Calculating sum of received pulses
                p_q = GaussianPulse((n_idx-1)/fs - tau_q, Vg, sig);
                gq_mn = alpha_q * p_q .* cos(2*pi*fc*((n_idx-1)/fs - tau_q));
                gq_mn = awgn(gq_mn, SNR_dB, 'measured');
                gq_prime(m_idx, n_idx) = gq_mn;
            end
        end
    end
    m = 1:M; n = 1:N; q = 1:Q;

    % Carrier signal
    [~,s] = ComplexCarrier(Vg, N, fc, fs);

    % Determining y^H(m,n) as the low-pass filtered version of p_q(mn)*s(n)
    yH_mn_unf = gq_prime .* s;
    yH_mn = lowpass(gq_prime .* s, fp, fs);
end

% Fast Fourier Transform function, with real and imaginary parts
% Amplitude not converted to mV
function [f_half, Y_mag_unf, Y_mag] = FFT(Tp, M, N, fc, fs, Vg, sig, c, fp, SNR_dB, range_bin, Tw)
    % Get signal at specific range bin
    [~, ~, ~, ~, yH_mn, yH_mn_unf] = ADC_MultiTarget(M, N, Tp, fc, fp, fs, Vg, sig, c, SNR_dB);
    yH_m_unf = yH_mn_unf(:, range_bin);  % Slow-time signal (M x 1)
    yH_m = yH_mn(:, range_bin); % Slow-time signal (M x 1)

    % Convert window duration to number of pulses
    W = round(Tw / Tp);
    if mod(W, 2) == 0
        W = W + 1;  % Make it odd for symmetric window
    end

    % Extract centered window in slow time
    m_center = round(M / 2);
    half_W = floor(W / 2);
    idx_start = max(1, m_center - half_W);
    idx_end = min(M, m_center + half_W);
    yH_window_unf = yH_m_unf(idx_start:idx_end);
    yH_window = yH_m(idx_start:idx_end);

    % Remove DC component
    yH_window_unf = yH_window_unf - mean(yH_window_unf);
    yH_window = yH_window - mean(yH_window);

    % Apply Hamming window
    yH_window_unf = yH_window_unf .* hamming(length(yH_window_unf));
    yH_window = yH_window .* hamming(length(yH_window));

    % FFT
    M_fft = 2^nextpow2(4 * W);
    f = linspace(0, 1/Tp, M_fft);  % Slow-time frequency axis (Hz)
    yH_pad_unf = [yH_window_unf', zeros(1, M_fft - length(yH_window_unf))];
    yH_pad = [yH_window', zeros(1, M_fft - length(yH_window))];
    YH_unf = fft(yH_pad_unf);
    YH = fft(yH_pad);
    Y_mag_unf = abs(YH_unf(1:M_fft/2)); % Magnitude spectrum (single-sided)
    Y_mag = abs(YH(1:M_fft/2));
    f_half = f(1:M_fft/2); % Corresponding frequency axis
end

% PMR calculation for multiple range bins
% Notice the result is out of scale since the FFT's ampltidue is not in mV
function [range, pmr] = PMR(Tp, M, N, fc, fs, Vg, sig, c, fp, SNR_dB, Tw)
    % FFT
    Q = 50;
    range = linspace(0.6, 4, Q);
    pmr = zeros(1, Q);
    for i = 1:Q
        [~, ~, Y_mag] = FFT(Tp, M, N, fc, fs, Vg, sig, c, fp, SNR_dB, i, Tw);
        pmr(i) = max(Y_mag) / median(Y_mag);
    end
end

% CWT Scalogram, with real and imaginary parts, windowed to slow time
function [wt, mag_wt, freq] = CWT(Tp, M, N, fc, fp, fs, Vg, sig, c, SNR_dB, range_bin, Tw)
    % Get signal at specific range bin (filtered)
    [~, ~, ~, ~, yH_mn, ~] = ADC_MultiTarget(M, N, Tp, fc, fp, fs, Vg, sig, c, SNR_dB);
    yH_m = yH_mn(:, range_bin);  % Slow-time signal

    % Window length in pulses
    W = round(Tw / Tp);
    if mod(W, 2) == 0
        W = W + 1;
    end

    % Extract centered window
    m_center = round(M / 2);
    half_W = floor(W / 2);
    idx_start = max(1, m_center - half_W);
    idx_end = min(M, m_center + half_W);
    yH_window = yH_m(idx_start:idx_end);

    % Remove DC and apply Hamming window
    yH_window = yH_window - mean(yH_window);
    yH_window = yH_window .* hamming(length(yH_window));

    % CWT
    fb = cwtfilterbank('SignalLength', length(yH_window), 'SamplingFrequency', 2.8, ...
    'FrequencyLimits', [0.05 5], 'Wavelet','amor'); % Changing the frequency limits due to function limitation

    [wt, freq] = cwt(yH_window, 'FilterBank', fb); % Changing the length, so that peaks at zero don't interfere

    mag_wt = sqrt(wt(:,:,1).^2 + wt(:,:,2).^2);
end

% MAD Score, with real and imaginary parts, windowed to slow time
function [mad_score, t_window] = MAD(Tp, M, N, fc, fp, fs, Vg, sig, c, SNR_dB, range_bin, Tw)
    % Get signal at specific range bin (filtered)
    [~, ~, ~, ~, yH_mn, ~] = ADC_MultiTarget(M, N, Tp, fc, fp, fs, Vg, sig, c, SNR_dB);
    yH_m = yH_mn(:, range_bin);  % Slow-time signal (complex)

    % Window length in pulses
    W = round(Tw / Tp);
    if mod(W, 2) == 0
        W = W + 1;
    end

    % Extract centered window
    m_center = round(M / 2);
    half_W = floor(W / 2);
    idx_start = max(1, m_center - half_W);
    idx_end = min(M, m_center + half_W);
    yH_window = yH_m(idx_start:idx_end);

    % Remove DC and apply Hamming window
    yH_window = yH_window - mean(yH_window);
    yH_window = yH_window .* hamming(length(yH_window));
    
    t_window = (0:length(yH_window)-1) * Tp;

    % Compute magnitude
    mag_yH = abs(yH_window);

    % Compute MAD in a sliding window
    winSize = 4;  % number of slow-time pulses in MAD window (can adjust)
    mad_score = movmad(mag_yH, winSize, 1);
end


% Variables ===============================================================
% Values should be later adjusted to make sense in the physical world

[N, M] = deal(623, 600); % Number of pulses, Total number of samples in each pulse interval
fs = 10000;
fc = 100; % Carrier frequency (in Hz)
Tp = 0.1; % Pulse repetition interval (PRI, in sec)
t = 1:1/fs:M*Tp; % Time span, arbitrary (in sec)
fast_time = (0:N-1)/fs * 1e6; % (in msec)
slow_time = (0:M-1) * Tp; % (in sec)
Vg = 1; % Amplitude
sig = 0.01; % Standard deviation
c = physconst('LightSpeed'); % Speed of ligth (in m/s)
R = 1.2; % Distance to the target (in meters)
fp = fs/100; % Passband frequency, arbitrary == only produces relevant results for some specific values
SNR_dB = 25; % Signal-to-noise ratio (in dB)
range_bins = [10, 50, 100, 150, 200, 250]; % Arbitrary
r = 100; % Specific range bin to perform FFT
Tw = 30; % Time window to perform FFT (in sec)
thold = 1000; % Arbitrary threshold for PMR selection
freqs = linspace(0.07, 10, length(slow_time)); % List of wavelet bases from min to max, same amount of points as time signal

% Main ====================================================================
% Single pulse
p = GaussianPulse(t, Vg, sig);

% Complex carrier
[nc, s] = ComplexCarrier(Vg, N, fc, fs);

% Multi-point model for the received signal
[m, n, q, gq_prime, yH_mn, yH_mn_unf] = ADC_MultiTarget(M, N, Tp, fc, fp, fs, Vg, sig, c, SNR_dB);

% FFT of multi-point model
[f_half, Y_mag_unf, Y_mag] = FFT(Tp, M, N, fc, fs, Vg, sig, c, fp, SNR_dB, r, Tw);

% PMR values for the plot
[range, pmr] = PMR(Tp, M, N, fc, fs, Vg, sig, c, fp, SNR_dB, Tw);

% CWT Scalogram of waveform sequence, regardless of PMR value (for now)
[wt, mag_wt, f] = CWT(Tp, M, N, fc, fp, fs, Vg, sig, c, SNR_dB, r, Tw);

% MAD Score
[mad_score, t_window] = MAD(Tp, M, N, fc, fp, fs, Vg, sig, c, SNR_dB, r, Tw);

% Plots ===================================================================

% Unfiltered |y^H(m,n)| (fast time on x, slow time on y)
figure;
imagesc(n, m, abs(yH_mn_unf));
xlabel('Fast Time Index (n)');
ylabel('Slow Time Index (m)');
title('Unfiltered |y^{H}(m,n)|');
%colormap hot;
colorbar;

% Filtered |y^H(m,n)| (fast time on x, slow time on y)
figure;
imagesc(n, m, abs(yH_mn));
xlabel('Fast Time Index (n)');
ylabel('Slow Time Index (m)');
title('Filtered |y^{H}(m,n)|');
%colormap hot;
colorbar;

% Slow time plot of the unffiltered signal
figure; hold on;
colors = lines(length(range_bins));
for i = 1:length(range_bins)
    n_idx = range_bins(i);
    yH_m_unf = yH_mn_unf(:, n_idx);
    plot(slow_time, yH_m_unf, 'Color', colors(i,:), 'DisplayName', ['n = ' num2str(n_idx)]);
end
xlabel('Slow time (t [s])');
ylabel('y^{H}(m, n)');
title('Examples of unfiltered signals');
legend show;
grid on;

% Slow time plot of the filtered signal
figure; hold on;
colors = lines(length(range_bins));
for i = 1:length(range_bins)
    n_idx = range_bins(i);
    yH_m = yH_mn(:, n_idx);
    plot(slow_time, yH_m, 'Color', colors(i,:), 'DisplayName', ['n = ' num2str(n_idx)]);
end
xlabel('Slow time (t [s])');
ylabel('y^{H}(m, n)');
title('Examples of filtered signals');
legend show;
grid on;

% FFT of both filtered and unfiltered slow time signals (for selected bin)
figure; hold on;
plot(f_half, Y_mag_unf, f_half, Y_mag);
xlabel('Frequency [Hz]');
ylabel('|FFT(y_{eq})|');
[Y_max, idx_max] = max(Y_mag);
maxstr = sprintf('Frequency of maximum for filtered signal = %.2f Hz', idx_max);
legend ('Unfiltered', 'Filtered', maxstr);
grid on;

% PMR for different ranges of the filtered frequencies
figure;
plot(range, pmr, '-o');
yline(thold);
xlabel('Range [m]');
ylabel('PMR');
title('PMR value to identify motion range bins');
legend('Values', 'Threshold');
grid on;
 
% CWT of signal above threshold (for selected bin, r)
figure;
imagesc(slow_time, f, abs(mag_wt));
xlabel('Time [s]');
ylabel('Frequency [Hz]');
title('CWT Scalogram for r = 100');
% colormap hot;
colorbar;

% CWT of signal over time (for selected bin, r)
figure; 
plot(slow_time, abs(mag_wt(5,:)));
xlabel('Time [s]');
ylabel('|CWT|');
grid on;

% MAD score of signal over time (for selected bin, r)
figure;
plot(t_window, mad_score);
xlabel('Time [s]');
ylabel('MAD Score');
title('MAD Score for r = 100');
grid on;