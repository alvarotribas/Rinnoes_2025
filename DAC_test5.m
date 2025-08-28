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
    space = abs(-0.99*sin(4*R-8));
    space(space < 0.25) = 0; % Half-wave retification, just an artificial trick to obtain the result for now
    
    % Attenuation model in free-space
    attenuation = 1 ./ (R.^1.7);
    idx = find(R >= 1, 1, 'first'); % Finding index where attenuation will begin
    attenuation(1:idx-1) = 1;

    % Behavior of the gain in the room
    alpha = space .* attenuation;
end

% Multi-point model
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
fp = fs/100; % Passband frequency, arbitrary == only produces relevant results for some specific values
SNR_dB = 35; % Signal-to-noise ratio (in dB)
range_bins = [10, 50, 100, 150, 200, 250]; % Arbitrary
up_factor = 10; % DAC parameter
Q = N; % Number of points in the model
R = linspace(0.6, 4.25, Q); % Distribution of distances

% Functions ===============================================================

% p(t)
p = GaussianPulse(t, Vg, sig);

% alpha(R)
[alpha, space, attenuation] = Gain(R);

% gs_prime(m,n)
gs_prime = ADC(Tp, M, fc, Vg, sig, c, fs, N, R, Q, alpha);

% Plots ===================================================================

% alpha(R)
figure;
plot(R, space, R, alpha, R, attenuation);
xlabel('R'); ylabel('\alpha'); title('Gain profile in the model');
legend('Space propagation', 'Attenuation', 'Result');
grid on;

% gs_prime(m,n,R)
figure;
imagesc(1:N, 1:M, abs(gs_prime));
xlabel('Fast time index (n)');
ylabel('Slow time index (m)');
title('ADC |g_{s}^{p}(m,n)|');
colorbar;

