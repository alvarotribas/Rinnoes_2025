% Notes
% - This code implements the solution done in DAC_test4.m into the code of
% slprof_test2.m, that is, the multi-point scenario
% - Besides the combination, it also addresses the range bin selection,
% through a formula that relates it to distance and sampling frequency

% Functions ===============================================================

% Generate a single Gaussian Pulse (p(t))
function p = GaussianPulse(t, Vg, sig)
    % Signal
    p = Vg*exp(-(t.^2) / (2 * sig^2));
end

% Gain profile in the model
function [alpha, space, attenuation] = Gain(R)
    % Model of the propagtion in the room (five peaks)
    %space = abs(-0.99*sin(4*R-8));
    %space(space < 0.25) = 0; % Half-wave retification, just an artificial trick to obtain the result for now
    width = 0.72;
    centers = linspace(0.75,3.9,5);
    space = zeros(size(R));
    for i = 1:length(centers)
        idx = (R >= centers(i)-width/2) & (R <= centers(i)+width/2);
        space(idx) = 0.99 * sin(pi * (R(idx) - (centers(i)-width/2)) / width);
    end

    % Attenuation model in free-space
    attenuation = 1 ./ (R.^2);
    idx = find(R >= 1, 1, 'first'); % Finding index where attenuation will begin
    attenuation(1:idx-1) = 1;

    % Behavior of the gain in the room
    alpha = space .* attenuation;
end

% Complex carrier signal
function [n,s] = ComplexCarrier(Vg, N, fc, fs)
    n = 1:N;
    s = Vg * exp(-1i*2*pi*fc*(n/fs));
end

% Signal received with multi-point contribution
function [t_rel, g_prime] = Received(t, Tp, M, fc, Vg, sig, alpha, c, R, Q)

    % Parameters
    g_prime = zeros(size(t));
    Ts = t(2) - t(1); % assume uniform sampling
    N = round(Tp / Ts); % samples per pulse

    % Respiratory movement model
    chest = 5; 
    a_all = normpdf(linspace(1, 5, Q), chest, 1.5); % Movement constant, normal distributed around torax region
    b_all = normpdf(linspace(0.1, 1, Q), chest/10, 0.2); % Movement frequency, also normal around torax
    
    for q = 1:Q
        % Parameters of the received pulse
        R_q = R(q);
        alpha_q = alpha(q);
        a = a_all(q);
        b = b_all(q);

        for m_idx = 1:M
            % Time window for single pulse
            t_start = (m_idx-1)*Tp;
            t_window = t - t_start;
    
            % Center of the fast-time window
            t_center = (N-1)/2 * Ts;
    
            % Target delay at pulse center
            t_abs_center = t_start + t_center;
            x_q = a*sin(b*t_abs_center); % motion of the chest
            tau_center = 2*(R_q + x_q)/c;
    
            % Fast-time relative to pulse center
            t_rel = t_window - t_center - tau_center;
    
            % Gaussian pulse
            p_rf = GaussianPulse(t_rel, Vg, sig);
    
            % Carrier modulation at relative fast-time
            g_prime = g_prime + alpha_q * p_rf .* cos(2*pi*fc*t_rel);
        end
    end
end

% Multi-point model
function y_mn = ADC(Tp, M, fc, Vg, sig, c, fs, N, R, Q, alpha)

    % Matrix with dimensions pulses x samples
    gs_single = zeros(M, N);
    gs_prime = zeros(M, N);

    % Respiratory movement model
    chest = 5; 
    a_all = normpdf(linspace(1, 10, Q), chest, 2); % Movement constant, normal distributed around torax region
    b_all = normpdf(linspace(0.1, 1, Q), chest/10, 0.2); % Movement frequency, also normal around torax

    for q = 1:Q
        % Parameters of the received pulse
        R_q = R(q);
        alpha_q = alpha(q);
        a = a_all(q);
        b = b_all(q);

        for m_idx = 1:M
            % Pulse start time
            t_start = (m_idx-1)*Tp;
    
            % Target motion (example sinusoidal)
            x_center = a*sin(b*t_start);
            tau_center = 2*(R_q + x_center)/c; % round-trip delay
    
            for n_idx = 1:N
                % Absolute fast-time
                t_fast = (n_idx-1)/fs;
    
                % Relative time wrt the delay
                t_rel = t_fast - tau_center;
    
                % Gaussian pulse
                p_mn = GaussianPulse(t_rel, Vg, sig);
    
                % Carrier modulation at relative fast-time
                gs_single(m_idx, n_idx) = alpha_q * p_mn * cos(2*pi*fc*t_rel);
            end

        end

        % Marking the correct range bin
        r = round(2*R_q*fs/c);
        gs_prime(:, r) = gs_single(:, r);

    end

    % Carrier signal
    [~,s] = ComplexCarrier(Vg, N, fc, fs);

    % Signal after demodulation
    y_mn = gs_prime .* s;

end

% DAC conversion
function [t_dac, g_dac] = DAC(y_mn, fs, up_factor, Tp)

    % Flatten pulses to chronological 1-D sequence
    g1d = y_mn.'; % N x M
    g1d = g1d(:); % (M*N) x 1

    % Resample once using sinc-like interpolation
    g_dac = resample(g1d, up_factor, 1); 

    % Time vector for reconstructed waveform
    fs_up = fs * up_factor;
    t_dac = (0:length(g_dac)-1).' / fs_up;

    % Correction to ADC signal through zero-padding
    g_dac(end+1) = 0;
    t_dac(end+1) = 0;

    % Attempt to filter the signal
    % bpFilt = designfilt('bandpassfir', 'FilterOrder', 450, 'CutoffFrequency1', 6.75e6, ...
    % 'CutoffFrequency2', 9.25e6, 'SampleRate', 1/Tp);

    g_dac = real(bandpass(g_dac, [7.75e6 8.25e6], 1/Tp));
    % g_dac = filtfilt(bpFilt, real(g_dac));
end

% Fast Fourier Transform function, with amplitude
function [f_half, Gprime_mag, Gdac_mag] = FFT(g_prime, g_dac, Tp, Tw, M)

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
    gprime_window = g_prime(idx_start:idx_end);
    gdac_window = g_dac(idx_start:idx_end);

    % Remove DC component
    gprime_window = gprime_window - mean(gprime_window);
    gdac_window = gdac_window - mean(gdac_window);

    % Apply Hamming window
    gprime_window = gprime_window .* hamming(length(gprime_window));
    gdac_window = gdac_window .* hamming(length(gdac_window));

    % FFT
    M_fft = 2^nextpow2(4 * W);
    f = linspace(0, 1/Tp, M_fft);  % Slow-time frequency axis (Hz)
    gprime_pad = [gprime_window', zeros(1, M_fft - length(gprime_window))];
    gdac_pad = [gdac_window', zeros(1, M_fft - length(gdac_window))];
    % gdac_pad = bandpass(gdac_pad, [7.5e6 8.5e6], 1/Tp);
    Gprime = fft(gprime_pad);
    Gdac = fft(gdac_pad);
    Gprime_mag = abs(Gprime(1:M_fft/2));
    Gdac_mag = abs(Gdac(1:M_fft/2));
    f_half = f(1:M_fft/2); % Corresponding frequency axis
end

% Comparing g_dac to g_prime via Hilbert and Fourier transforms
function [env1, env2, G1, G2, f1, f2] = comp_Hilbert_FFT(fs1, g1, fs2, g2)
    % Alignment of the signals
    L = min(length(g1), length(g2));
    g1 = g1(1:L);
    g2 = g2(1:L);

    % Hilbert transform + normalization
    env1 = abs(hilbert(g1)); env1 = env1 / max(env1); 
    env2 = abs(hilbert(g2)); env2 = env2 / max(env2);

    % FFT (positive frequencies only)
    N1 = length(g1); N2 = length(g2);
    G1_full = fft(g1, N1);
    G2_full = fft(g2, N2);
    
    % Indices for positive frequencies [0, fs/2)
    idx1 = 1:floor(N1/2);
    idx2 = 1:floor(N2/2);
    
    G1 = abs(G1_full(idx1));
    G2 = abs(G2_full(idx2));
    
    f1 = (0:floor(N1/2)-1) * (fs1/N1);
    f2 = (0:floor(N2/2)-1) * (fs2/N2);

end

% Variables ===============================================================

[N, M] = deal(623, 600); % Number of pulses, Total number of samples in each pulse interval
fs = 16e9; % Sampling frequency (in Hz)
Ts = 1/fs; % Sampling period (in sec)
fc = 5e9; % Carrier frequency (in Hz)
Tp = N/fs; % Pulse repetition interval (PRI, in sec)
t = 0:1/fs:M*Tp; % Time span, arbitrary
fast_time = (0:N-1)/fs * 1e6; % (in msec)
slow_time = (0:M-1) * Tp; % (in sec)
Vg = 1; % Amplitude
sig = 1e-9; % Standard deviation
c = physconst('LightSpeed'); % Speed of ligth (in m/s)
fp = fs/5; % Passband frequency, arbitrary == only produces relevant results for some specific values
SNR_dB = 35; % Signal-to-noise ratio (in dB)
range_bins = [10, 50, 100, 150, 200, 250]; % Arbitrary
up_factor = 1; % DAC parameter
Q = N; % Number of points in the model
R = linspace(0.6, 4.25, Q); % Distribution of distances
fs1 = 1/(t(2)-t(1)); fs2 = fs*up_factor; % Frequencies for the comparative FFT
Tw = 30*1e-6; % Time window to perform FFT 

% Functions ===============================================================

% p(t)
p = GaussianPulse(t, Vg, sig);

% alpha(R)
[alpha, space, attenuation] = Gain(R);

% g_prime
[t_rel, g_prime] = Received(t, Tp, M, fc, Vg, sig, alpha, c, R, Q);

% gs_prime(m,n)
y_mn = ADC(Tp, M, fc, Vg, sig, c, fs, N, R, Q, alpha);

% r(t)
[t_dac, g_dac] = DAC(y_mn, fs, up_factor, Tp);

% Manual FFT
[f_half, Gprime_mag, Gdac_mag] = FFT(g_prime', g_dac, Tp, Tw, M);

% H[g(t)], F[g(t)]
[env1, env2, G1, G2, f1, f2] = comp_Hilbert_FFT(fs1, g_prime, fs2, g_dac);

% Plots ===================================================================

% % alpha(R)
% figure;
% plot(R, space, 'b', R, attenuation, 'g', R, alpha, 'r');
% xlabel('R'); ylabel('\alpha'); title('Gain profile in the model');
% legend('Space propagation', 'Attenuation', 'Result');
% grid on;

% g_prime(t)
figure;
plot(t*1e6, g_prime);
xlabel('Time [\mus]'); ylabel('Amplitude');
title('Received signal g^p(t)');
grid on;

% gs_prime(m,n,R)
figure;
imagesc(1:N, 1:M, abs(y_mn));
xlabel('Fast time index (n)');
ylabel('Slow time index (m)');
title('ADC |g_{s}^{p}(m,n)| with demodulation');
colorbar;

% DAC
figure;
plot(t_dac*1e6, g_dac);
xlabel('Time [\mus]'); ylabel('Amplitude');
title('r(t) reconstruction after the DAC');
grid on;

% Manual FFT
figure;
plot(f_half, 20*log10(Gprime_mag/max(Gprime_mag)), 'b', f_half, 20*log10(Gdac_mag/max(Gdac_mag)), 'r');
xlabel('Frequency [Hz]');
ylabel('20 \log_{10}(G_{mag}/\max(G_{mag}))');
legend('Received g\_prime', 'Reconstructed g\_dac');
grid on;

% H[g(t)], F[g(t)]
figure; subplot(2,1,1);
plot(t*1e6, env1, 'b', t_dac*1e6, env2, 'r');
xlabel('Time [\mus]');
ylabel('Envelope amplitude');
legend('Received g\_prime', 'Reconstructed g\_dac');
title('Envelope comparison');
grid on;

subplot(2,1,2);
plot(f1/1e6, 20*log10(G1/max(G1)), 'b', f2/1e6, 20*log10(G2/max(G2)), 'r');
xlabel('Frequency [MHz]');
ylabel('Magnitude [dB]');
legend('Received g\_prime', 'Reconstructed g\_dac');
title('Spectral comparison');
grid on;

% Testing =================================================================