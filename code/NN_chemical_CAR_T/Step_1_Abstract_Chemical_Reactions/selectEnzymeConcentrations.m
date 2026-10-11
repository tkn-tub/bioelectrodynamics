function D = selectEnzymeConcentrations(K, dt, WX, opt)
%SELECTENZYMECONCENTRATIONS  Algorithm 1 of the paper: enzyme concentrations
%   of a chemical differentiator from its timing and gain targets.
%
%   D = selectEnzymeConcentrations(E_Table, dt, WX)
%   D = selectEnzymeConcentrations(E_Table, dt, WX, opt)
%
%   E_Table  table from readKineticsTable with rows r_S, r_Sr, r_F, r_Fr, r_I
%            and columns KM [mM] and kcat [1/s]; or a struct with fields
%            K.KM.S, K.KM.Sr, K.KM.F, K.KM.Fr, K.KM.I and K.kcat.(same).
%   dt       relative delay between slow and fast branch, Delta t [s]
%   WX       bandwidth of the input signal [rad/s]
%   opt      optional struct, any of:
%     wS, wF   branch gains                                   (default 1, 1)
%     pF       folding probability                            (default 1)
%     Emax     largest feasible enzyme concentration [mM]     (default 0.01 = 10 uM)
%     kI       rate of r_I acting on the fast pool [1/s]      (default 0)
%     Xmax     largest input concentration [mM]; enables the linear-regime check
%     tauS     starting tau_S [s]                             (default 2*dt)
%     epsW     what "WX*tau_S << 1" means                     (default 0.1)
%     epsLin   what "<< K_M" means                            (default 0.1)
%     names    row names in E_Table  (default {'r_S','r_Sr','r_F','r_Fr','r_I'})
%
%   D has the concentrations E.S, E.Sr, E.F, E.Fr [mM], the timing tauS, tauF,
%   alpha, dt, the rate constants k.S, k.Sr, k.F, k.Fr [1/s], and D.feasible.
%
%   Every rate constant follows Eq. (12), k = pF*[E]*kcat/KM, so that
%   tau = 1/k_r and w = k_f/k_r (Eq. (2)).

    if nargin < 4, opt = struct(); end
    wS     = getOpt(opt, 'wS', 1);
    wF     = getOpt(opt, 'wF', 1);
    pF     = getOpt(opt, 'pF', 1);
    Emax   = getOpt(opt, 'Emax', 0.01);
    kI     = getOpt(opt, 'kI', 0);
    Xmax   = getOpt(opt, 'Xmax', []);
    epsW   = getOpt(opt, 'epsW', 0.1);
    epsLin = getOpt(opt, 'epsLin', 0.1);
    names  = getOpt(opt, 'names', {'r_S', 'r_Sr', 'r_F', 'r_Fr', 'r_I'});

    % ---- Input: kcat_i, KM_i for i in {S, Sr, F, Fr, I} --------------------
    if ~isstruct(K)
        tag = {'S', 'Sr', 'F', 'Fr', 'I'};
        T = K;  K = struct();
        for i = 1:numel(tag)
            K.KM.(tag{i})   = T{names{i}, 'KM'};
            K.kcat.(tag{i}) = T{names{i}, 'kcat'};
        end
    end
    KM = K.KM;  kcat = K.kcat;
    v = [KM.S KM.Sr KM.F KM.Fr kcat.S kcat.Sr kcat.F kcat.Fr];
    assert(all(isfinite(v)) && all(v > 0), 'KM and kcat of r_S, r_Sr, r_F and r_Fr must be positive numbers.');
    assert(dt > 0 && WX > 0, 'dt and WX must be positive.');

    % ---- 1: choose tau_S with WX*tau_S << 1 and tau_S > dt -----------------
    tauSmax = epsW/WX;                       % largest tau_S that still follows the input
    tauS    = getOpt(opt, 'tauS', 2*dt);
    assert(tauS > dt, 'tau_S must be larger than dt.');
    if dt >= tauSmax
        warning('selectEnzymeConcentrations:bandwidth', ...
            'dt = %.3g s is not below %.3g s, the largest tau_S that satisfies WX*tau_S <= %.2g.', dt, tauSmax, epsW);
    end

    % Linear regime (line 8, second condition). It depends only on the input
    % level and the gains, not on tau_S or dt, so it is checked once.
    linear = true;
    if ~isempty(Xmax)
        need = [Xmax, Xmax, wS*Xmax, wF*Xmax];          % substrate seen by r_S, r_F, r_Sr, r_Fr
        lim  = epsLin*[KM.S, KM.F, KM.Sr, KM.Fr];
        who  = {'r_S', 'r_F', 'r_Sr', 'r_Fr'};
        if isfinite(KM.I), need(end+1) = wF*Xmax; lim(end+1) = epsLin*KM.I; who{end+1} = 'r_I'; end
        bad = need > lim;
        if any(bad)
            linear = false;
            warning('selectEnzymeConcentrations:linear', ...
                ['Not in the linear regime for %s. Reduce the input level below %.3g mM ' ...
                 'or the gains; changing tau_S or dt does not help.'], ...
                strjoin(who(bad), ', '), min(lim./(need/Xmax)));
        end
    end

    % ---- 2-11: repeat ... until all concentrations are feasible ------------
    feasible = false;
    for it = 1:1000
        alpha = 1 - dt/tauS;                                  % 3: tau_F = alpha*tau_S
        tauF  = alpha*tauS;
        E.Sr  = KM.Sr/(pF*kcat.Sr*tauS);                      % 4: slow-branch removal
        kFr   = 1/tauF - kI;                                  % 5: fast pool removed by r_Fr and r_I
        if kFr <= 0
            kFr = 0;                                          % r_I alone is already fast enough
        end
        E.Fr  = KM.Fr*kFr/(pF*kcat.Fr);                       % 6: fast-branch removal
        E.S   = wS*KM.S/(pF*kcat.S*tauS);                     % 7: k_f = w*k_r
        E.F   = wF*KM.F/(pF*kcat.F*tauF);
        if all([E.S E.Sr E.F E.Fr] <= Emax)                   % 8
            feasible = true;  break
        end
        % 9: increase tau_S while the input is still followed, then reduce dt
        if tauS < tauSmax
            tauS = min(1.1*tauS, tauSmax);
        elseif dt > 1e-3*tauS
            dt = 0.9*dt;
        else
            break
        end
    end
    if kFr == 0 && kI > 0
        warning('selectEnzymeConcentrations:kI', ...
            'k_I = %.3g 1/s already exceeds 1/tau_F = %.3g 1/s, so [E]_Fr = 0 and the fast branch is faster than requested.', kI, 1/tauF);
    end
    if ~feasible
        warning('selectEnzymeConcentrations:infeasible', ...
            'No feasible set below Emax = %.3g mM; returning the closest one (tau_S = %.3g s, dt = %.3g s).', Emax, tauS, dt);
    end

    D.E = E;
    D.tauS = tauS;  D.tauF = tauF;  D.alpha = alpha;  D.dt = dt;
    D.k = struct('S', pF*E.S*kcat.S/KM.S, 'Sr', pF*E.Sr*kcat.Sr/KM.Sr, ...
                 'F', pF*E.F*kcat.F/KM.F, 'Fr', pF*E.Fr*kcat.Fr/KM.Fr, 'I', kI);
    D.wS = D.k.S/D.k.Sr;
    D.wF = D.k.F/(D.k.Fr + kI);
    D.WXtauS = WX*tauS;
    D.feasible = feasible && linear;

    fprintf('Algorithm 1: tau_S = %.3g s, tau_F = %.3g s, dt = %.3g s, alpha = %.3g, WX*tau_S = %.3g\n', ...
        tauS, tauF, dt, alpha, WX*tauS);
    fprintf('  [E]_S  = %9.4g uM   k_S  = %.4g 1/s\n', 1e3*E.S,  D.k.S);
    fprintf('  [E]_Sr = %9.4g uM   k_Sr = %.4g 1/s   (gain w_S = %.3g)\n', 1e3*E.Sr, D.k.Sr, D.wS);
    fprintf('  [E]_F  = %9.4g uM   k_F  = %.4g 1/s\n', 1e3*E.F,  D.k.F);
    fprintf('  [E]_Fr = %9.4g uM   k_Fr = %.4g 1/s   (gain w_F = %.3g, k_I = %.3g 1/s)\n', 1e3*E.Fr, D.k.Fr, D.wF, kI);
    fprintf('  feasible below %.4g uM: %d;  linear regime: %d\n', 1e3*Emax, feasible, linear);
end

function v = getOpt(opt, name, default)
    v = default;
    if isfield(opt, name) && ~isempty(opt.(name)), v = opt.(name); end
end
