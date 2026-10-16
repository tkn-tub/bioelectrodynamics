function [m_d1, aux] = D1_IFFL(m_in, dt, p)
% D1_IFFL  Block-1 differentiator (inverted polarity) via an incoherent
%          feed-forward loop. m_d1 rests near zero and pulses up when m_in falls.
%
%   [m_d1, aux] = D1_IFFL(m_in, dt)        uses default parameters
%   [m_d1, aux] = D1_IFFL(m_in, dt, p)     override any field of p (see below)
%
% Topology (reaction nodes r_*, species on the arrows):
%   fast arm : m_in --r_Ff--> m_Ff --r_Fd--> m_Fd --r_Fr--> 0     (repressor)
%   slow arm : m_in --r_Sf--> m_Sf --r_Sd--> m_d1                  (produces output)
%   r_I      : m_d1 --(enzyme = m_Fd, SATURATED in m_d1)--> 0      (incoherent repression)
%   r_lk     : m_d1 --> 0                                          (small leak)
%
% The subtraction is done by r_I: the fast-delayed species m_Fd acts as the
% enzyme that degrades the output m_d1, run in its saturated (zero-order-in-m_d1)
% regime so it removes a FLUX set by the fast arm without adding to the output's
% effective rate constant. Hence low baseline AND preserved dynamics, unlike the
% bimolecular annihilator.

    m_in = m_in(:);  N = numel(m_in);  m0 = m_in(1);

    % ---------- parameters (override via p) ----------
    if nargin < 3, p = struct(); end
    def = struct( ...
        'tau_F',30, 'tau_S',60, 'tau_r',30, ...   % fast/slow delays, feed/removal time [s]
        'Vfac',1.5, ...                            % r_I max rate / resting flux  (>1 pulls baseline down)
        'K_I',2e-6, ...                            % r_I Michaelis const [mM]  (sets the output floor)
        'k_leak',1/120, ...                        % output leak [1/s]
        'burnin_s',3000);                          % settle time to the resting state [s]
    fn = fieldnames(def);
    for i = 1:numel(fn), if ~isfield(p,fn{i}), p.(fn{i}) = def.(fn{i}); end, end

    kFf = 1/p.tau_r;  kSf = kFf;      % equal feed into both arms
    kFd = 1/p.tau_F;  kSd = 1/p.tau_S;
    kFr = 1/p.tau_r;                  % fast repressor removal
    flux0 = kSf*m0;                   % resting production of the output
    V_I   = p.Vfac*flux0;             % saturated max rate of r_I

    % ---------- resting pools ----------
    mFf = kFf*m0/kFd;  mSf = kSf*m0/kSd;  mFd = kFd*mFf/kFr;  mFd_ref = mFd;
    Z   = p.K_I;

    % ---------- burn-in to the true resting state ----------
    for k = 1:ceil(p.burnin_s/dt)
        rep = V_I*(mFd/mFd_ref)*Z/(p.K_I+Z);
        mFf = mFf + dt*(kFf*m0 - kFd*mFf);
        mSf = mSf + dt*(kSf*m0 - kSd*mSf);
        mFd = mFd + dt*(kFd*mFf - kFr*mFd);
        Z   = max(Z + dt*(kSd*mSf - rep - p.k_leak*Z), 0);
    end

    % ---------- run ----------
    m_d1 = zeros(N,1);  m_Fd = zeros(N,1);
    for n = 1:N
        m_d1(n) = Z;  m_Fd(n) = mFd;
        rep = V_I*(mFd/mFd_ref)*Z/(p.K_I+Z);            % r_I, saturated in m_d1
        mFf = mFf + dt*(kFf*m_in(n) - kFd*mFf);
        mSf = mSf + dt*(kSf*m_in(n) - kSd*mSf);
        mFd = mFd + dt*(kFd*mFf     - kFr*mFd);
        Z   = max(Z + dt*(kSd*mSf - rep - p.k_leak*Z), 0);
    end

    aux = struct('m_Fd',m_Fd,'baseline',m_d1(1),'V_I',V_I, ...
                 'E_rI_uM',1e3*V_I/0.909125);          % LDH k_cat; feasibility check vs E_max=20 uM
end
