%% Feature_extraction_NN_no_rDp.m
%
% Differentiator simplified to a 3-reaction IFFL (no r_Dp):
%   r_DF: emptyset -> m_d            (fast tracker IS the output)
%   r_DS: emptyset -> m_Ds           (slow tracker)
%   r_Di: m_d + m_Ds -> sink         (standard mass-action; both consumed)
%
% This matches the COBRA-mapping topology decision (m_Ds stoichiometrically
% consumed in r_Di, S(m_Ds, r_Di) < 0).
%
% NOTES vs previous version:
%   - r_Dp removed; m_Df renamed to m_d.
%   - FLAG-1 dissolved: with r_Dp gone, the "m_Df*m_Ds vs m_d*m_Ds" ambiguity
%     disappears (they were the same species; now they actually are).
%   - m_Ds is now consumed by r_Di (stoichiometric). If ps is too small the
%     pool depletes; bump ps if you see the bump fail after one event.
%   - FLAG-3 still applies: m_d and m_Ds have no first-order decay.
%     m_d returns to zero via r_Di only. m_Ds returns to zero via r_Di only.
%     For multi-pulse use, add small kf_deg/kS_deg or rely on r_Di to clean up.

% clc; clear; close all;

fontsize = 15;

%% =========================
% Parameters
%==========================

% --- Differentiator block (D) ---
% r_DF: emptyset -> m_d  (fast tracker, also the output)
pf      = 1;
k_DF    = 2;
K_DF    = 0.5;
k_d_deg = 0;     % first-order decay of m_d (kept zero; r_Di does the cleanup)

% r_DS: emptyset -> m_Ds  (slow tracker)
ps      = 0.99;
k_DS    = 2;
K_DS    = 0.5;
k_Ds_deg = 0;    % first-order decay of m_Ds (kept zero)

% r_Di: m_d + m_Ds -> sink  (mass-action, both substrates consumed)
k_Di    = 20;

% --- Counter block (C) ---
V_Cp    = 1.0;   K_Cp    = 0.05;
k_Cr    = 1.0;
k_Cg    = 2.0;
m_Cs_amp = 0.009;

% --- NN block (N) ---
V_Ni = 1.0;  K_Ni = 0.1;  n_Ni = 2;
V_Na = 1.0;  K_Na = 0.1;  n_Na = 2;
V_Nb = 1.0;  K_Nb = 0.1;  n_Nb = 2;
V_No = 1.0;  K_No = 0.1;  n_No = 2;
m_Nc_const = 1.0;
k_Ni_deg = 1.0; k_Na_deg = 1.0; k_Nb_deg = 1.0; k_out_deg = 1.0;

% --- Communication channel ---
M_0  = 1e4;
D_co = 66e-13;
d    = 20e-6;
Tp   = 2;
Ts   = 0.01;
fs   = 1/Ts;
sequence = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
r_rx = 0.15e-6;
V_rx = 4/3*pi*r_rx^3;

%% =========================
% Time vector and Tx waveform
%==========================
samples_per_pulse = floor(Tp * fs);
total_samples     = length(sequence) * samples_per_pulse;
t   = (0:total_samples-1)' / fs;
mTx = repelem(sequence, samples_per_pulse)' * M_0;

%% =========================
% Diffusion channel -> received concentration m_in
%==========================
t_h = t; t_h(1) = 1/fs;
h   = V_rx ./ (4*pi*D_co*t_h).^(3/2) .* exp(-d^2 ./ (4*D_co*t_h));
m_in = conv(h, mTx);
m_in = m_in(1:length(t));

%% =========================
% State arrays
%==========================
N = length(t);

% Differentiator block
m_d  = zeros(N,1);    % fast tracker AND output (cross-block signal)
m_Ds = zeros(N,1);    % slow tracker

% Counter block
m_c  = zeros(N,1);
m_Cs = m_Cs_amp * ones(N,1);

% NN block
m_Ni = zeros(N,1);
m_Na = zeros(N,1);
m_Nb = zeros(N,1);
m_out = zeros(N,1);

% Diagnostic
v_DF_arr = zeros(N,1);
v_DS_arr = zeros(N,1);
v_Di_arr = zeros(N,1);

%% =========================
% Discrete-time simulation
%==========================
for n = 1:N-1

    dt = t(n+1) - t(n);

    % ============================================================
    % DIFFERENTIATOR BLOCK   (3-reaction IFFL, no r_Dp)
    % ============================================================
    % r_DF: emptyset -> m_d
    v_DF = pf * k_DF * ( m_in(n) / (K_DF + m_in(n)) );
    v_DF_arr(n+1) = v_DF;

    % r_DS: emptyset -> m_Ds
    v_DS = ps * k_DS * ( m_in(n) / (K_DS + m_in(n)) );
    v_DS_arr(n+1) = v_DS;

    % r_Di: m_d + m_Ds -> sink   (mass-action; consumes both)
    v_Di = k_Di * m_d(n) * m_Ds(n);
    v_Di_arr(n+1) = v_Di;

    % Stoichiometric updates (each species written once, contributions summed)
    m_d(n+1)  = max( m_d(n)  + dt * ( v_DF - v_Di - k_d_deg  * m_d(n) ),  0 );
    m_Ds(n+1) = max( m_Ds(n) + dt * ( v_DS - v_Di - k_Ds_deg * m_Ds(n) ), 0 );

    % ============================================================
    % COUNTER BLOCK
    % ============================================================
    % r_Cp:  m_d -> m_c   (catalytic in m_d; m_d not subtracted here)
    v_Cp = V_Cp * m_d(n) / (K_Cp + m_d(n));
    v_Cr = k_Cr * m_c(n) * m_Cs(n);
    v_Cg = k_Cg * m_d(n) * m_Cs(n);

    m_c(n+1)  = max( m_c(n)  + dt * ( v_Cp - v_Cr ), 0 );
    m_Cs(n+1) = max( m_Cs(n) - dt*v_Cg + m_Cs_amp,   0 );

    % ============================================================
    % NN BLOCK
    % ============================================================
    act_Ni = (m_c(n)^n_Ni) / (K_Ni^n_Ni + m_c(n)^n_Ni);
    v_Ni   = V_Ni * act_Ni;
    m_Ni(n+1) = max( m_Ni(n) + dt*( v_Ni - k_Ni_deg * m_Ni(n) ), 0 );

    act_Na = (m_Ni(n)^n_Na) / (K_Na^n_Na + m_Ni(n)^n_Na);
    v_Na   = V_Na * act_Na;
    m_Na(n+1) = max( m_Na(n) + dt*( v_Na - k_Na_deg * m_Na(n) ), 0 );

    act_Nb = (m_Ni(n)^n_Nb) / (K_Nb^n_Nb + m_Ni(n)^n_Nb);
    v_Nb   = V_Nb * act_Nb;
    m_Nb(n+1) = max( m_Nb(n) + dt*( v_Nb - k_Nb_deg * m_Nb(n) ), 0 );

    drive_No = m_Nc_const * (m_Na(n) + m_Nb(n));
    act_No   = (drive_No^n_No) / (K_No^n_No + drive_No^n_No);
    v_No     = V_No * act_No;
    m_out(n+1) = max( m_out(n) + dt*( v_No - k_out_deg * m_out(n) ), 0 );

end

%% =========================
% Plot 1: m_in and m_d
%==========================
figure;
plot(t, m_in, 'LineWidth', 2); hold on;
ylabel('$m_\mathrm{in}$ concentration', 'Interpreter', 'latex');
yyaxis right
plot(t, m_d, 'LineWidth', 2);
ylabel('$m_d$ concentration', 'Interpreter', 'latex');
xlabel('Time [s]', 'Interpreter', 'latex');
legend({'$m_\mathrm{in}$', '$m_d$'}, 'Interpreter', 'latex', 'Location', 'best');
set(gca, 'FontSize', fontsize); grid on;
title('Differentiator (3-rxn IFFL): $m_\mathrm{in}\rightarrow m_d$', 'Interpreter', 'latex');

%% =========================
% Plot 2: differentiator internals (debug)
%==========================
figure;
plot(t, m_d,  'LineWidth', 2); hold on;
plot(t, m_Ds, 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Concentration', 'Interpreter', 'latex');
legend({'$m_d$ (fast)', '$m_{Ds}$ (slow)'}, 'Interpreter', 'latex', 'Location', 'best');
set(gca, 'FontSize', fontsize); grid on;
title('IFFL internals', 'Interpreter', 'latex');

%% =========================
% Plot 3: Counter
%==========================
figure;
plot(t, m_in, 'LineWidth', 2);
ylabel('$m_\mathrm{in}$', 'Interpreter', 'latex');
yyaxis right
plot(t, m_c, 'LineWidth', 2);
ylabel('$m_c$', 'Interpreter', 'latex');
xlabel('Time [s]', 'Interpreter', 'latex');
legend({'$m_\mathrm{in}$', '$m_c$'}, 'Interpreter', 'latex', 'Location', 'best');
set(gca, 'FontSize', fontsize); grid on;
title('Counter: $m_d\rightarrow m_c$', 'Interpreter', 'latex');

%% =========================
% Plot 4: NN response
%==========================
figure;
plot(t, m_c,  'LineWidth', 2); hold on;
plot(t, m_Ni, 'LineWidth', 2);
plot(t, m_Na, 'LineWidth', 2);
plot(t, m_Nb, 'LineWidth', 2);
plot(t, m_out,'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Concentration', 'Interpreter', 'latex');
legend({'$m_c$', '$m_{Ni}$', '$m_{Na}$', '$m_{Nb}$', '$m_\mathrm{out}$'}, ...
       'Interpreter', 'latex', 'Location', 'best');
set(gca, 'FontSize', fontsize); grid on;
title('NN: $m_c\rightarrow m_\mathrm{out}$', 'Interpreter', 'latex');
