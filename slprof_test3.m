% Notes:
% - This code is written to test the sleep apnea detection algorithm,
% following the results from a paper. No specficic range bin has been
% selected according to the criteria stablished in slprof_test2.m.
% - A later code will be written combining the actual results of all parts.

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

% Correlation function
function detection = SleepApneaDetection(M, N, y_mn, epsilon)
    % Definition
    X = M/50; % Number of consecutive samples, arbitrary
    L = 10/(X*10^(-6)); % Length of apnea duration interval
    
    % Windowing
    chi_q = zeros(size(N));
    for k = 1:N
        yi = y_mn(:,k);
        chi_M = max(yi(1:X));
        chi_q = chi_q + chi_M;
    end

    rho_bar = zeros(size(N));
    for k = L+1:N
        q = k-L;
        rho_M = max(chi_q(q:k));
        rho_m = min(chi_q(q:k));
        rho_bar = rho_bar + (1/2)*(rho_M + rho_m);
    end

    % Sleep apnea alarm flags
    alarm = zeros(size(length(rho_bar)));
    for i = 1:length(rho_bar)
        if ((1-epsilon)*rho_bar(i) <= chi_q(i)) && (chi_q(i) <= (1+epsilon)*rho_bar(i))
            alarm = alarm + 1;
        end
    end

    % Sleep apnea detection algorithm
    target_value = 1; mask = ones(1,L); % Flag number and mask to identify it
    logical_op = alarm == target_value;

    detection = conv(double(logical_op), mask, 'valid') == L;
end

% Variables ===============================================================
[N, M] = deal(623, 10000); % Number of pulses, Total number of samples in each pulse interval
fs = 10000;
fc = 100; % Carrier frequency (in Hz)
Tp = 0.1; % Pulse repetition interval (PRI, in sec)
t = 1:1/fs:M*Tp; % Time span, arbitrary (in sec)
fast_time = (0:N-1)/fs * 1e6; % (in msec)
slow_time = (0:M-1) * Tp/100; % (in sec)
Vg = 1; % Amplitude
sig = 0.01; % Standard deviation
c = physconst('LightSpeed'); % Speed of ligth (in m/s)
R = 1.2; % Distance to the target (in meters)
fp = fs/100; % Passband frequency, arbitrary == only produces relevant results for some specific values
SNR_dB = 20; % Signal-to-noise ratio (in dB)
range_bins = [10, 50, 100, 150, 200, 250]; % Arbitrary
r = 100; % Specific range bin to perform FFT
Tw = 30; % Time window to perform FFT (in sec)
thold = 1000; % Arbitrary threshold for PMR selection
freqs = linspace(0.07, 10, length(slow_time)); % List of wavelet bases from min to max, same amount of points as time signal
alpha = 0.9;
epsilon = 0.1;

% Main ====================================================================
% Received signals at Rx (g^p)
g_prime = Received(t, Tp, M, fc, Vg, sig, alpha, c, R);

% Sleep apnea detection algorithm
detection = SleepApneaDetection(M, N, g_prime, epsilon);
if any(detection)
    disp('Interval of the apnea ='); disp(find(detection));
else
    fprintf('No apnea was detected');
end

% Plots ===================================================================
% Plot of examples of range bins in the filtered baseband signal
figure;
plot(t, g_prime);
xlabel('t [s]');
ylabel('Amplitude');
legend('Real');
grid on;

