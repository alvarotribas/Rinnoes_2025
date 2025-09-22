% Notes
% - This code is a more robust test of the sleep apnea detection algorithm,
% where I attempt to use the results in DAC_test4.m -- particularly the DAC
% signal over time -- to find the correct alarms. Previously, I only
% guaranteed the algorithm works, but didn't have correct values for the
% signal parameters.
% - To do this, i attempt to make the selection algorithm based on arbitrary
% criteria, to check if the reconstruction of the signal holds. In the final
% version, this will be done with the FFT and MAD criteria stablished in
% slprof_test2.m.

% Functions ===============================================================

% Generate a single Gaussian Pulse (p(t))
function p = GaussianPulse(t, Vg, sig)
    % Signal
    p = Vg*exp(-(t.^2) / (2 * sig^2));
end

% Received RF pulse signal
function [t_rel, g_prime] = Received(t, Tp, M, fc, Vg, sig, alpha, c, R, Q)

    % Parameters
    g_prime = zeros(size(t));
    Ts = t(2) - t(1); % assume uniform sampling
    N = round(Tp / Ts); % samples per pulse

    for q = 1:Q

        R_q = R(q);
        alpha_q = alpha(q);

        for m_idx = 1:M
            % Time window for single pulse
            t_start = (m_idx-1)*Tp;
            t_window = t - t_start;
    
            % Center of the fast-time window
            t_center = (N-1)/2 * Ts;
    
            % Target delay at pulse center
            t_abs_center = t_start + t_center;
            x_center = 2*sin(1*t_abs_center); % motion of the chest
            tau_center = 2*(R_q + x_center)/c;
    
            % Fast-time relative to pulse center
            t_rel = t_window - t_center - tau_center;
    
            % Gaussian pulse
            p_rf = GaussianPulse(t_rel, Vg, sig);
    
            % Carrier modulation at relative fast-time
            g_prime = g_prime + alpha_q * p_rf .* cos(2*pi*fc*t_rel);
        end

    end

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

% Multi-point model, ADC conversion
function gs_prime = ADC(Tp, M, fc, Vg, sig, c, fs, N, R, Q, alpha)

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

% Correlation function
function [detection, chi_q, samples, alarm, apnea_intervals] = SleepApneaDetection(y_i, epsilon)
    
    % Parameters
    l = length(y_i); % Length of the received signal
    X = round(l/100000); % Number of consecutive samples, arbitrary
    L = 1000; % Length of apnea duration interval, in number of samples. Arbitrary
    samples = floor(l/X); % Number of X samples in the signal
    

    % Windowing
    chi_q = zeros(1, samples);
    for j = 1:samples
        a = (j-1)*X + 1; 
        b = min(j*X, l); 
        chi_q(j) = max(y_i(a:b));
    end

    % Calculating the average correlation
    rho_bar = zeros(1, samples);
    for j = L+1:samples
        k = j-L;
        rho_M = max(chi_q(k:j));
        rho_m = min(chi_q(k:j));
        rho_bar(j) = (1/2)*(rho_M + rho_m);
    end

    % Sleep apnea flag condition
    alarm = zeros(1, samples);
    for i = 1:samples
        if rho_bar(i) > 0
            if ((1-epsilon)*rho_bar(i) <= chi_q(i)) && (chi_q(i) <= (1+epsilon)*rho_bar(i))
                alarm(i) = 1;
            end
        end
    end

    % Sleep apnea detection algorithm
    mask = ones(1,L);
    streaks = conv(double(alarm), mask, 'valid') >= L;

    detection = zeros(1, samples);
    detection(L:end) = streaks;

    % Storing the indexes of the apnea intervals
    delta = diff([0 detection 0]);
    start_idx = find(delta==1);
    end_idx = find(delta==-1) -1;

    apnea_intervals = [(start_idx-1)*X + 1, end_idx*X]; 
    apnea_intervals(apnea_intervals > l) = l; % clip to signal length
end

% Variables ===============================================================

[N, M] = deal(623, 600); % Number of pulses, Total number of samples in each pulse interval
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
up_factor = 10; % DAC parameter
Q = N; % Number of points in the model
R = linspace(0.6, 4.25, Q); % Distribution of distances
epsilon = 0.7;

% Main ====================================================================

%g_prime(t)
[t_rel, g_prime] = Received(t, Tp, M, fc, Vg, sig, alpha, c, R, Q);

% alpha(R)
[alpha, space, attenuation] = Gain(R);

% gs_prime(m,n)
gs_prime = ADC(Tp, M, fc, Vg, sig, c, fs, N, R, Q, alpha);

% Reconstruct waveform
[t_dac, g_dac, g_dac_matrix] = DAC(gs_prime, fs, up_factor);

% Sleep apnea detection algorithm
[detection, chi_q, samples, alarm, apnea_intervals] = SleepApneaDetection(g_dac, epsilon);

if isempty(apnea_intervals)
    disp('No apnea detected');
else
    disp('Apnea intervals (start sample, end sample):');
    disp(apnea_intervals);
end

% Plots ===================================================================

% g_prime(t)
figure;
plot(t*1e3, g_prime);
xlabel('Time [ms]'); ylabel('Amplitude');
title('Received signal g^p(t)');

% gs_prime(m,n,R)
figure;
imagesc(1:N, 1:M, abs(gs_prime));
xlabel('Fast time index (n)');
ylabel('Slow time index (m)');
title('ADC |g_{s}^{p}(m,n)|');
colorbar;

% DAC
figure;
plot(t_dac*1e6, g_dac);
xlabel('Time [\mus]'); ylabel('Amplitude');
title('Samples of pulses reconstructed by the DAC');

% Alarm plot, to see where apnea was detected
figure;
plot(1:samples, chi_q); hold on;
plot(1:samples, alarm);
xlabel('Samples');
ylabel('Amplitude');
legend('Correlation', 'Alarm');
grid on;