% Notes:
% - Correct the code such that the signals are in adequate time frames
% - For simplicity, no range bin is being selected. The single point model
% is being adopted as if it were a mutiple point model, with all bins.

% Functions ===============================================================
% Generate a single Gaussian Pulse (p(t))
function p = GaussianPulse(t, Vg, sig)
    % Signal
    p = Vg*exp(-(t.^2) / (2 * sig^2));
end

% Generate the baseband signal
function g = Baseband(t, Tp, M, fc, Vg, sig)
    % Repeated pulses
    g = zeros(size(t));

    for m = 1:M
        tau_m = t - m*Tp;
        p_bb = GaussianPulse(tau_m, Vg, sig);
        g = g + p_bb .* cos(2*pi*fc*tau_m);
    end
end

% Received RF pulse signal
function g_prime = Received(t, Tp, M, fc, Vg, sig, alpha, c, R)
    % Round-trip time
    x = 2*sin(1*t); % Movement of the chest, estimated
    tau = 2*(R + x)/c; % Round-trip time, is very small due to speed of light

    g_prime = zeros(size(t));

    for m = 1:M
        tau_m = t - m*Tp;
        p_rf = GaussianPulse(tau_m - tau, Vg, sig);
        g_prime = g_prime + alpha * p_rf .* cos(2*pi*fc*(tau_m-tau));
    end
end

% Analog-to-digital converter (ADC)
function [n, m, gs_prime] = ADC(Tp, M, fc, Vg, sig, alpha, c, R, fs, N, SNR_dB)
    % Inspired by the function above
    Ts = 1/fs;
    gs_prime = zeros(M,N);
    
    for m_idx = 1:M
        for n_idx = 1:N
            t_mn = (m_idx-1)*Tp + (n_idx-1)*Ts;
            x_mn = 1*sin(1*t_mn); % Estimated, could be changed
            tau_mn = 2*(R+x_mn)/c;

            p_mn = GaussianPulse((n_idx-1)/fs - tau_mn, Vg, sig);
            gs_mn = alpha * p_mn .* cos(2*pi*fc*((n_idx-1)/fs - tau_mn));
            gs_mn = awgn(gs_mn, SNR_dB, 'measured');

            % Multivariable equation
            gs_prime(m_idx,n_idx) = gs_mn;
        end
    end
    m = 1:M;
    n = 1:N;
end

% Complex carrier signal
function [n,s] = ComplexCarrier(Vg, N, fc, fs)
    n = 1:N;
    s = Vg * exp(-1i*2*pi*fc*(n/fs));
end

% Complex complete baseband signal, in the scenario with a single point
function [n, m, y_mn] = LPF(Tp, M, N, fc, fs, Vg, sig, alpha, c, R, fp, SNR_dB)
    % Load the equations that are multiplied 
    [~, ~, gs_prime] = ADC(Tp, M, fc, Vg, sig, alpha, c, R, fs, N, SNR_dB);
    [~, s] = ComplexCarrier(Vg, N, fc, fs);

    % Matrix multiplication
    y_mn = gs_prime .* s;

    % Filtering
    y_mn = lowpass(y_mn, fp, fs);
    n = 1:N;
    m = 1:M;
end

% Complex baseband signal from the equation (to compare with LPF)
function [n, m, y_eq] = LPFEquation(Tp, M, N, fc, fs, Vg, sig, alpha, c, R)
    Ts = 1/fs;
    y_eq = zeros(M,N);

    for m = 1:M
        for n = 1:N
            t_mn = m*Tp + n*Ts;
            x_mn = 0.5*sin(5*t_mn); % Estimated, could be changed
            tau_mn = 2*(R+x_mn)/c;
            %Equation
            y_eq_aux = (alpha*Vg/2)*exp(-(n*Ts-tau_mn)^2/(2*sig^2))*exp(-1i*2*pi*fc*tau_mn);
            
            % Building the heatmap from the equation
            y_eq(m,n) = y_eq_aux;
        end
    end

    % Defining variables for the future
    m = 1:M;
    n = 1:N;
end

% Fast Fourier Transform function, with real and imaginary parts
function [f_half, Y_mag] = FFT(Tp, M, N, fc, fs, Vg, sig, alpha, c, R, fp, SNR_dB, range_bin, Tw)
    % Get signal at specific range bin
    [~, ~, y_mn] = LPF(Tp, M, N, fc, fs, Vg, sig, alpha, c, R, fp, SNR_dB);
    y_m = y_mn(:, range_bin);  % Slow-time signal (M x 1)

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
    y_windowed = y_m(idx_start:idx_end);

    % Remove DC component
    y_windowed = y_windowed - mean(y_windowed);

    % Apply Hamming window
    y_windowed = y_windowed .* hamming(length(y_windowed));

    % FFT
    M_fft = 2^nextpow2(4 * W);
    f = linspace(0, 1/Tp, M_fft);  % Slow-time frequency axis (Hz)
    y_padded = [y_windowed', zeros(1, M_fft - length(y_windowed))];
    Y_full = fft(y_padded);
    Y_mag = abs(Y_full(1:M_fft/2)); % Magnitude spectrum (single-sided)
    f_half = f(1:M_fft/2); % Corresponding frequency axis
end

% Variables ===============================================================
% Values should be later adjusted to make sense in the physical world

[N, M] = deal(1000, 300); % Number of pulses, Total number of samples in each pulse interval
fs = 8e9; % Sampling frequency (in Hz)
fc = 4e9; % Carrier frequency (in Hz)
Tp = 4e-3; % Pulse repetition interval (PRI, in sec)
t = 1:1/fs:M*Tp; % Time span, arbitrary
fast_time = (0:N-1)/fs * 1e6; % (in msec)
slow_time = (0:M-1) * Tp; % (in sec)
Vg = 1; % Amplitude
sig = 9e-9; % Standard deviation
c = physconst('LightSpeed'); % Speed of ligth (in m/s)
R = 1.2; % Distance to the target
alpha = 0.95; % Channel gain
fp = fs/100; % Passband frequency, arbitrary == only produces relevant results for some specific values
SNR_dB = 50; % Signal-to-noise ratio (in dB)
range_bins = [10, 50, 100, 150, 200, 250]; % Arbitrary

% Main ====================================================================
% Single pulse
p = GaussianPulse(t, Vg, sig);

% Baseband signal
g = Baseband(t, Tp, M, fc, Vg, sig);

% Received RF pulse signal
g_prime = awgn(Received(t, Tp, M, fc, Vg, sig, alpha, c, R), SNR_dB, 'measured');

% Signal after ADC
%[n, m, gs_prime] = ADC(Tp, M, fc, Vg, sig, alpha, c, R, fs, N, SNR_dB);

% Complex carrier
%[nc, s] = ComplexCarrier(Vg, N, fc, fs);

% Complex complete baseband from definition
%[n_def, m_def, y_mn] = LPF(Tp, M, N, fc, fs, Vg, sig, alpha, c, R, fp, SNR_dB);

% Complex complete baseband from equation
%[n_cbb, m_cbb, y_eq] = LPFEquation(Tp, M, N, fc, fs, Vg, sig, alpha, c, R);

% FFT of arbitrary range bin
%[f_half, Y_mn] = FFT(Tp, M, N, fc, fs, Vg, sig, alpha, c, R, fp, SNR_dB, range_bin);

% Saving the most important functions
save("g_prime_4GHz.mat", "g_prime", "t");
%save("adc_4GHz.mat", "gs_prime"); % Remember to plot the absolute value later
%save("lpf_4GHz.mat", "y_mn"); % Remeber to plot the absolute value later

% % Plots ===================================================================
% % Original train of pulses
% figure;
% plot(t, p);
% xlabel('Time [s]');
% ylabel('p(t)');
% grid on;
% 
% % Received signal
% figure;
% plot(t, g_prime);
% xlabel('t [s]');
% ylabel('Amplitude');
% legend('Real');
% grid on;
% 
% % Plot of the received pulse (fast vs slow time)
% figure;
% imagesc(1:length(n), 1:length(m), abs(gs_prime));
% xlabel('Fast time index (n)');
% ylabel('Slow time index (m)');
% title('ADC |g_{s}^{p}(m,n)|');
% colorbar;
% 
% % Plot of the complex baseband signal, filtered (fast vs slow time)
% figure;
% imagesc(n_def, m_def, abs(y_mn));
% xlabel('Fast time index (n)');
% ylabel('Slow time index (m)');
% title('LPF |y_{mn}(m,n)|');
% colorbar;
% 
% % Plot of examples of range bins in the filtered baseband signal
% figure; hold on;
% colors = lines(length(range_bins));
% for i = 1:length(range_bins)
%     n_idx = range_bins(i);
%     y_m = y_mn(:, n_idx);
%     plot(slow_time, y_m, 'Color', colors(i,:), 'DisplayName', ['n = ' num2str(n_idx)]);
% end
% xlabel('t [s]');
% ylabel('y_{mn}(m, n)');
% title('Examples of filtered range bins');
% legend show;
% grid on;
% 
% % Plot of FFT for a specific range bin
% figure; hold on;
% bins = 100;  
% colors2 = lines(length(bins));
% for i = 1:length(bins)
%     elem = bins(i);
%     [f_half, Y_mn] = FFT(Tp, M, N, fc, fs, Vg, sig, alpha, c, R, fp, SNR_dB, elem, 30);
%     plot(f_half, Y_mn);
% end
% xlabel('Frequency [Hz]');
% % xlim([0 2]);
% ylabel('|FFT(y_{eq})|');
% legend ('Absolute');
% grid on;