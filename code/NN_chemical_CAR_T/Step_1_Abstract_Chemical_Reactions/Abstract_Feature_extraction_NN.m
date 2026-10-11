% BLOCK 0: selected adapter path 61 from Block_7_Selected_Motifs(4).txt.
% cAMP --PDE1--> AMP --NTD7--> adenosine --PUNP1--> ribose-1-P
%       --PYNP2r (reverse)--> uridine --r1156--> UMP = m_in.
% cAMP is a dynamic pool, produced from buffered ATP by an ADNCYC anchor.
% EV-dependent activation of ADNCYC is an ASSUMED one-site transduction law;
% the receptor/extracellular signaling steps are not resolved here.
% Vesicle counts are an activation signal, never a chemical substrate.
% Block 0 Km, kcat, enzyme amounts, cofactors and initial pools are illustrative.
% The report states KineticsValidated=0 and RegulatoryEvidenceComplete=0.
% This simulates a candidate UMP input; it does not infer a measured cellular input.
% Directional apparent MM rates follow the selected path; reverse backfluxes
% and feedback by other cellular reactions are not modeled.
% ATP, Pi, uracil and dUTP are held constant; H2O is implicit.
% UMP is NOT buffered: slow and fast SOURCE reactions withdraw it.
% D1 is updated to two lag branches and a shared canceller (five reactions).
% D2/counter laws are retained; NN now implements the affine sum of Eq. (6).
% Unrelated motif reactions in the report
% are not substituted for this working abstract design.
% Main for loop: named concentrations and explicit Euler; no pool solver.
% Initial pools are at equilibrium for the first EV sample to avoid a
% pulse caused solely by starting the adapter with empty intermediates.
% Original 13 figure roles retained; figures 14/15 expose Block 0 dynamics.
%
% Branch notation: each differentiator has local reactions r_S1, r_S2, r_F.
% D1 / D2 suffixes distinguish the stages in MATLAB variable names.
% D1 slow upper branch: X --r_S1--> M1 --r_S2--> Y1.
% D1 fast lower branch: X --r_F--> F1_delay --r_R--> products.
% D1 shared cancellation: F1_delay + Y1 --r_I--> products.
% m_S remains the SLOW memory; m_d1 remains the physical output Y1.
% m_F is the NEW fast memory. Existing signal names are retained.
% Both producer/consumer pairs are lag elements. r_I couples their pools,
% so their loaded relaxation times are not ideal independent time delays.
% t and every tau below are in SECONDS. Main output plots use t/60 [minutes].
% timing_profile='current' preserves the previous reaction calibration.
% timing_profile='faster' reduces memory lag while preserving UMP source fluxes.
% The faster profile also lowers pulse amplitude; inspect the counter response.
% Tau edits here recalibrate assumed pools/Km values, not measured affinities.
% D2 slow upper branch: Y1 --r_S1--> M2 --r_S2--> products.
% D2 fast lower branch: Y1 --r_F--> Y2; r_I consumes M2 + Y2.
% D2 r_R removes Y2 for extra output cleanup.
% Examples: V_S / v_S; V_S1_D2 / v_S1_D2.
% D1 capacities and equilibria are recalibrated for its five reactions.
% D2 and counter reaction laws are retained.
% NN: C and a maintained bias B produce the SAME pool Z=m_Ni.
% dZ/dt = k_Nw*C+k_Nb*B-k_Ni_deg*Z; Z_ss = w_N*C+b_N.
% The default identity readout m_out=Z implements the affine distance estimate.
% Optional Hill activation adds output production and removal after Z.
% Counter C regulates production without being consumed by the NN readout;
% buffered precursors are implicit, as in the previous NN source reactions.
% A literal C -> Z conversion would also need a sink in the counter balance.
% m_Na/m_Nb are contribution diagnostics for Z, not extra chemical pools.
%
% Feature extraction: falling-slope stage -> positive-pulse stage -> counter -> NN.
% Named concentration arrays and explicit Euler updates in the main for loop.
% Stage 1 now has FIVE reactions: r_S1, r_S2, r_F, r_R and r_I.
% Delayed slow production competes with the filtered fast canceller.
% The unchanged downstream D2 sources also consume m_d1 in consuming mode.
% Stage 2 now has FIVE reactions (four for positive adaptation, one cleanup):
%   r_F_D2: Y1 -> Y2                  direct output production
%   r_S1_D2: Y1 -> M2                  slow inhibitor production
%   r_S2_D2: M2 -> products            inhibitor turnover
%   r_I_D2: M2 + Y2 -> products       consuming inhibition/removal
%   r_R_D2: Y2 -> products            extra output cleanup
% Rising Y1 raises Y2 before M2 catches up. Falling Y1 leaves excess M2,
% which removes the output pulse. This implements positive polarity.
% r_R_D2 allows fast cleanup of a large pulse but makes adaptation approximate.
%
% IMPORTANT MODEL LIMITS:
% Continuous production leaves a finite basal Y2, not exact zero.
% The defaults below illustrate strong saturation: Km_I_D2_Y=1e-13 mM.
% This is NOT a measured enzyme constant or a validated metabolic mapping.
% The small Y1 excursion requires substantial amplification to make the
% pulse large compared with basal Y2. Verify affinities/mechanisms before
% interpreting this example as a physical enzyme implementation.
% Basal Y2=1e-9 mM is ~0.6 molecules in 1 pL: this deterministic example
% also needs concentration/volume scaling and stochastic assessment for cells.
% The four-reaction core is smaller, but alone can have a long output tail.
%
% Only Y2 activates the counter; ordinary one-site kinetics (n_gate2=1).
% Y2 is a nonconsuming regulator of buffered counter production in this
% abstraction. No direct Y1-counter connection. A mechanism must be mapped.
% Counter is a leaky analog accumulator, not a held elapsed-time register.
% Diagnostic slopes/rectification are used ONLY for comparison plots.
%
% Concentrations and Km: mM; V: mM/s; kcat: 1/s; time: s.
% Factorized two-substrate MM is an assumed apparent rate law.
% No implicit pool solver or concentration clipping is used.
% The fast second stage uses visible Euler substeps inside the for loop.
% Other blocks use one explicit Euler step per input sample.
% fzero is used ONLY to initialize the constant-input equilibrium.

clc;
clearvars;
close all;

% Including parent folders up to your project folder named 'code'.
addingPathParentFolderByName('code'); %[output:6f7ea38d]
Parameters
%%
%[text] #### Loading Enzymes' dataset
%[text] Description: This block uses your credentials at "Brenda\_credentials.txt" for searching at [Brenda website](https://www.brenda-enzymes.info/functional_parameter_statistics.php) the enzyme parameters. Your credentials should be included in a txt file with two lines, the first line includes your email, and the second line your password. The enzymes are identified from the file "Block\_7\_Selected\_Motifs.txt", which is completed by the second step and includes the chemical reaction network of the algorithm. After the code reads the enzyme parameters it prepares an excel table with the listed values per enzyme.
%Uncomment below for reading the constants from https://www.brenda-enzymes.info/functional_parameter_statistics.php 

% %Loggin in Brenda for quering the enzyme parameters https://www.brenda-enzymes.info/functional_parameter_statistics.php 
% % Default values
% credentialsFile = fullfile(pwd,'BRENDA_credentials.txt');
% if isfile(credentialsFile)
%     fid = fopen(credentialsFile,'rt','n','UTF-8');
%     assert(fid ~= -1, 'Cannot open credentials file.');
% 
%     emailFromFile    = fgetl(fid);
%     passwordFromFile = fgetl(fid);
%     fclose(fid);
% 
%     assert(ischar(emailFromFile) && ischar(passwordFromFile), ...
%         'Credentials file must contain two lines.');
%     assert(~isempty(strtrim(emailFromFile)) && ~isempty(passwordFromFile), ...
%         'Email and password must not be empty.');
% 
%     email    = strtrim(emailFromFile);
%     password = passwordFromFile; % Preserve password spaces.
% 
%     clear fid emailFromFile passwordFromFile
% end
% 
% 
% brendaKineticsTable('Block_7_Selected_Motifs.txt', 'kinetics_BRENDA.xlsx', email, password);

%Reading the table of enzyme kinetic parameters
E_Table = readKineticsTable('kinetics_BRENDA.xlsx', 1e-3);   % [E] = 1e-3 mM everywhere %[output:1f520a66] %[output:672d5d47] %[output:808e23bb] %[output:770b142c]

%%
%[text] ## Load sample data
load Dataset_cell2cell.mat
Dataset = rmmissing(Dataset);
dist = table2array(Dataset(:,1));
vesicles = table2array(Dataset(:,2:end));

row = 1;
gain_in = 1e-4;                     % OLD direct-input reference only [mM/vesicle]
vesicleInput = vesicles(row,:).';
t = (0:numel(vesicleInput)-1)' * dt0;
dtInput = dt0;
legacy_input_ref = gain_in * vesicleInput; % comparison only; never drives the chemistry
N = length(t);
assert(N >= 3 && all(isfinite(vesicleInput)) && all(vesicleInput >= 0) && vesicleInput(1) > 0, ...
    'Use a finite nonnegative EV signal with a positive first sample.');

%%
%[text] ## Block 0 design concentrations and buffered cosubstrates
%[text] Concentration targets \[mM\], not measurements. Capacity calibration follows below, after the differentiator's basal UMP consumption is known.
cAMP_ref = 0.01;
AMP_ref = 0.02; adn_ref = 0.02; r1p_ref = 0.02; uri_ref = 0.02;
UMP_ref = 0.05;
u_ref = UMP_ref;                    % first-stage design input, before initialization
K_EV_AC = 10*vesicleInput(1);       % activation constant in DATASET EV units, not mM
m_ATP_B0 = 2; m_Pi_B0 = 1; m_uracil_B0 = 0.05; m_dUTP_B0 = 0.05;
Km_ADNCYC_ATP = 0.2;
Km_PDE1_cAMP = 10*cAMP_ref;
Km_NTD7_AMP = 10*AMP_ref;
Km_PUNP1_adn = 10*adn_ref; Km_PUNP1_Pi = 0.1;
Km_PYNP2r_r1p = 10*r1p_ref; Km_PYNP2r_uracil = 0.005;
Km_r1156_uri = 10*uri_ref; Km_r1156_dUTP = 0.005;
% Effective turnover assumptions [1/s]; replace with verified enzyme values.
kcat_ADNCYC = 100; kcat_PDE1 = 100; kcat_NTD7 = 100;
kcat_PUNP1 = 100; kcat_PYNP2r = 100; kcat_r1156 = 100;
slope = diff(vesicles,1,2)/dtInput;   % vesicles/s, plotting reference
[peak,loc] = min(slope,[],2);
peakTime = (loc-1)*dtInput;
%%
%[text] ## Parameters: first differentiator and positive-pulse second stage
timing_profile = 'current';    % choose 'current' or 'faster'
plot_timing_checks = 1;        % additional slope/output comparison figure
tau_M1 = 30; % slow memory / output relaxation targets [s]
tau_Y1 = 3;       
tau_F1 = 3;                   % NEW fast-pool relaxation, Y1 held fixed [s]
source_ratio_1 = 2;            % fast/slow source ratio; leaves fast turnover flux
Y1_ref = 1e-2;                % design first-stage basal concentration [mM]
tau_M2 = 10;                  % inhibitor relaxation with Y2 held fixed [s]
tau_Y2 = 1;                   % added output cleanup time at low Y2 [s]
Y2_ref = 1e-9;                % design second-stage basal concentration [mM]
J_F_D2_ref = 1e-4;              % direct basal production flux [mM/s]
source_ratio_2 = 2;           % v_S1_D2/v_F_D2; must exceed one
K_gate2 = 1.2e-8;             % ordinary Y2 counter-activation constant [mM]
% n_gate2=1: fractions are written directly without Hill powers.

assert(strcmpi(timing_profile,'current') || strcmpi(timing_profile,'faster'), ...
    'Choose timing_profile = current or faster.');
if strcmpi(timing_profile,'faster')
    tau_M1 = 15;              % less first-stage memory lag
    tau_F1 = 1;               % quicker fast cancelling pool
    tau_M2 = 5;               % quicker second-stage inhibitor memory
    % Keep tau_Y1=3 and tau_Y2=1: preserve the basal source/cleanup fluxes.
end

% These are abstract calibration targets, not BRENDA measurements.
Km_D2_X = 10*Y1_ref;           % SHARED input Km of r_F_D2 and r_S1_D2
M2_ref = source_ratio_2*J_F_D2_ref*tau_M2*(10/11);
Km_D2_M = 10*M2_ref;           % SHARED memory Km of r_S2_D2 and r_I_D2
Km_I_D2_Y = 1e-13;              % very small assumed affinity; see header
Km_R_D2_Y = 1e-3;
V_R_D2_target = Km_R_D2_Y/tau_Y2;
J_R_D2_ref = V_R_D2_target*Y2_ref/(Km_R_D2_Y+Y2_ref);
J_I_D2_ref = J_F_D2_ref-J_R_D2_ref;
J_S2_D2_ref = (source_ratio_2-1)*J_F_D2_ref+J_R_D2_ref;
assert(J_I_D2_ref > 0 && source_ratio_2 > 1, 'Invalid stage-2 design fluxes.');
V_F_D2_target = J_F_D2_ref/(Y1_ref/(Km_D2_X+Y1_ref));
V_S1_D2_target = source_ratio_2*V_F_D2_target;
V_S2_D2_target = J_S2_D2_ref/(M2_ref/(Km_D2_M+M2_ref));
V_I_D2_target = J_I_D2_ref/ ...
    ((M2_ref/(Km_D2_M+M2_ref))*(Y2_ref/(Km_I_D2_Y+Y2_ref)));

% FIVE-REACTION D1: slow memory, fast memory and shared cancellation.
% At the design equilibrium:
% v_S = v_Sk = v_I_D1 + downstream Y1 withdrawals;
% v_F = v_Fk + v_I_D1.
consume_interstage=1;
gm1 = 10/11; gy1 = 10/11; gf1 = 10/11;
load_derivative = consume_interstage*(1+source_ratio_2)*J_F_D2_ref * ...
    (Km_D2_X/(Km_D2_X+Y1_ref))/Y1_ref;
clearance_derivative1 = 1/tau_Y1-load_derivative;
assert(clearance_derivative1 > 0 && source_ratio_1 > 1, ...
    'Check first-stage time constants and fast/slow source ratio.');
J_I_D1 = clearance_derivative1*Y1_ref/gy1;
J_S1_D1 = J_I_D1+consume_interstage*(1+source_ratio_2)*J_F_D2_ref;
J_F_D1 = source_ratio_1*J_S1_D1;
J_R_D1 = J_F_D1-J_I_D1;
M1_ref = J_S1_D1*tau_M1*gm1;
F1_delay_ref = J_F_D1*tau_F1*gf1;

%===========================================
% Enzyme constants
%===========================================

%Differentiator 1
%KM constant
KM_S = E_Table{"r_S","KM"};%10*u_ref;
KM_Sr = E_Table{"r_Sr","KM"};% = 10*M1_ref;
KM_F = E_Table{"r_F","KM"};%KM_S;
KM_Fr = E_Table{"r_Fr","KM"};%10*F1_delay_ref;
KM_I = E_Table{"r_I","KM"};%KM_Fr;      % shared fast-pool Km permits scalar initialization

% Specific activity
% SA [umol/(min*mg)] * MW [kDa] /60 -> kcat [1/s].
kcat_S = E_Table{"r_S","kcat"};% 10^2.5*48.901/60;
kcat_Sr = E_Table{"r_Sr","kcat"};% kcat_S;
kcat_F = E_Table{"r_F","kcat"};% 10^1.5*75/60;
kcat_Fr = E_Table{"r_Fr","kcat"};%kcat_F;
kcat_I = E_Table{"r_I","kcat"};%kcat_F;

%Enzyme concentration
opt.Emax = 0.02; %Departing from the largest feasible concentration [mM] = 20 uM
D = selectEnzymeConcentrations(E_Table, 1, 1/60, opt);   % dt = 1 s, WX = 1/60 rad/s %[output:3c07477a]

D.tauS, D.tauF, D.alpha %[output:4bab2dfa] %[output:2358b501] %[output:0c8c3c4e]

V_S = kcat_S*D.E.S; %E_S1_D1;
V_Sr = kcat_Sr*D.E.Sr;%E_S2_D1;
V_F = kcat_F*D.E.F;%E_F_D1*1.5;
V_Fr = kcat_Fr*D.E.Fr;%E_R_D1;
% V_I = kcat_I*E_I_D1;
E_I = V_Fr/kcat_I;
V_I = kcat_I*E_I;% The inhibitor reaction is setup to have the same speed as the fast branch


 V_F_D2 = kcat_F_D2*E_F_D2; V_S1_D2 = kcat_S1_D2*E_S1_D2; %[output:2e69004f]
% V_S2_D2 = kcat_S2_D2*E_S2_D2; V_I_D2 = kcat_I_D2*E_I_D2; V_R_D2 = kcat_R_D2*E_R_D2;



% --- BRENDA-Km edit (1/4): design references that follow from the measured Km ---
% The Km values are now fixed by BRENDA, so the reference concentrations are
% derived from them (before: Km = 10 x reference).
% Input UMP: a tenth of the smaller Km keeps both UMP sources in their linear
% range, so their ratio stays constant when the input changes.
u_lin = min(KM_S,KM_F)/10;
% Relaxation rate of the UMP pool at basal flux; it is advanced with one
% explicit Euler step per sample, which needs dt0*rate <= 0.5.
UMP_pool_rate = @(u) J_S1_D1*KM_S./(u.*(KM_S+u)) + J_F_D1*KM_F./(u.*(KM_F+u));
u_ref = u_lin;
if dt0*UMP_pool_rate(u_lin) > 0.5
    u_ref = fzero(@(u) dt0*UMP_pool_rate(u)-0.5,[u_lin 1e3*max(KM_S,KM_F)]);
    warning(['u_ref raised from %.3g to %.3g mM to keep the one-step Euler update of UMP stable ' ...
        'at dt0 = %.3g s. r_S then works at %.0f%% saturation (r_F at %.0f%%).'], ...
        u_lin,u_ref,dt0,100*u_ref/(KM_S+u_ref),100*u_ref/(KM_F+u_ref));
end
% Slow memory M1: same basal flux J_S1_D1 and same relaxation time tau_M1 as
% before, now solved for the measured KM_Sr from
%   V*M/(KM+M) = J   and   V*KM/(KM+M)^2 = 1/tau_M1.
% For KM_Sr = 10*M1_ref this reduces to the old M1_ref = J_S1_D1*tau_M1*10/11.
M1_ref = (-KM_Sr+sqrt(KM_Sr^2+4*J_S1_D1*tau_M1*KM_Sr))/2;
% --- end of edit (1/4) ---

%KM_I = 10*Y1_ref;
input_fraction_ref = u_ref/(KM_S+u_ref);
% --- BRENDA-Km edit (5/7): r_Fr and r_I now have DIFFERENT fast-pool Km values ---
% Fast pool F1: same basal fluxes and same relaxation time tau_F1 as before,
% now solved for the two measured Km values from
%   J_R_D1*KM_Fr/(F*(KM_Fr+F)) + J_I_D1*KM_I/(F*(KM_I+F)) = 1/tau_F1.
% For KM_I = KM_Fr = 10*F1_delay_ref this reduces to the old J_F_D1*tau_F1*10/11.
F1_relaxation = @(F) J_R_D1*KM_Fr./(F.*(KM_Fr+F))+J_I_D1*KM_I./(F.*(KM_I+F))-1/tau_F1;
F1_delay_ref = fzero(F1_relaxation,[1e-12*KM_Fr 1e6*max(KM_Fr,KM_I)]);
fast_fraction_ref = F1_delay_ref/(KM_Fr+F1_delay_ref);
% --- BRENDA-Km edit (2/4): each UMP source is calibrated with its OWN Km ---
V_S1_D1_target = J_S1_D1/input_fraction_ref;            % r_S uses KM_S (unchanged)
V_S2_D1_target = J_S1_D1/(M1_ref/(KM_Sr+M1_ref));
V_F_D1_target = J_F_D1/(u_ref/(KM_F+u_ref));            % was J_F_D1/input_fraction_ref (assumed KM_F = KM_S)
% --- end of edit (2/4) ---
V_R_D1_target = J_R_D1/fast_fraction_ref;
% r_I is calibrated with its OWN Km (was fast_fraction_ref, i.e. KM_I = KM_Fr assumed).
V_I_D1_target = J_I_D1/ ...
    ((F1_delay_ref/(KM_I+F1_delay_ref))*(Y1_ref/(KM_I+Y1_ref)));
% --- end of edit (5/7) ---


kcat_F_D2 = kcat_S; kcat_S1_D2 = kcat_S;
kcat_S2_D2 = kcat_S; kcat_I_D2 = kcat_F; kcat_R_D2 = kcat_F;
E_S1_D1 = V_S1_D1_target/kcat_S; E_S2_D1 = V_S2_D1_target/kcat_Sr;
E_F_D1 = V_F_D1_target/kcat_F;
E_R_D1 = V_R_D1_target/kcat_Fr; E_I_D1 = V_I_D1_target/kcat_I;
E_F_D2 = V_F_D2_target/kcat_F_D2; E_S1_D2 = V_S1_D2_target/kcat_S1_D2;
E_S2_D2 = V_S2_D2_target/kcat_S2_D2; E_I_D2 = V_I_D2_target/kcat_I_D2;
E_R_D2 = V_R_D2_target/kcat_R_D2;
% assert(abs(V_S1_D2-source_ratio_2*V_F_D2) <= 1e-12*V_S1_D2, ...
%     'The positive-pulse calibration requires matched input source capacities.');
% reaction_parameters = table( ...
%     {'r_S1_D1';'r_S2_D1';'r_F_D1';'r_R_D1';'r_I_D1';'r_F_D2';'r_S1_D2';'r_S2_D2';'r_I_D2';'r_R_D2'}, ...
%     [KM_S;KM_Sr;KM_F;KM_Fr;KM_I;Km_D2_X;Km_D2_X;Km_D2_M;Km_D2_M;Km_R_D2_Y], ...
%     [NaN;NaN;NaN;NaN;KM_I;NaN;NaN;NaN;Km_I_D2_Y;NaN], ...
%     [kcat_S;kcat_Sr;kcat_F;kcat_Fr;kcat_I;kcat_F_D2;kcat_S1_D2;kcat_S2_D2;kcat_I_D2;kcat_R_D2], ...
%     [E_S1_D1;E_S2_D1;E_F_D1;E_R_D1;E_I_D1;E_F_D2;E_S1_D2;E_S2_D2;E_I_D2;E_R_D2], ...
%     [V_S;V_Sr;V_F;V_Fr;V_I;V_F_D2;V_S1_D2;V_S2_D2;V_I_D2;V_R_D2], ...
%     'VariableNames',{'Reaction','Km_first_mM','Km_second_mM', ...
%                      'kcat_per_s','Enzyme_mM','Vmax_mM_per_s'});
% disp(reaction_parameters);
%%
%[text] ## Block 0 capacities: match the basal UMP demand of the first stage
%[text] UMP is consumed by BOTH r\_S1\_D1 and r\_F\_D1 (1:1 in this model).
% --- Algorithm-1 edit: size the UMP supply for the linear range of r_S and r_F ---
% Before: J_B0_ref = J_S1_D1+J_F_D1;   (flux targets of the old calibration; with the
% Algorithm 1 capacities V_S and V_F it left m_in near 0.04 mM, i.e. r_S saturated)
% Now the supply is the withdrawal of r_S and r_F at the target input level,
% a tenth of the smaller Km, so that both reactions work in their linear range.
m_in_target = min(KM_S,KM_F)/10;
J_B0_ref = V_S*m_in_target/(KM_S+m_in_target)+V_F*m_in_target/(KM_F+m_in_target);
% --- end of Algorithm-1 edit ---
EV_activity_base = vesicleInput(1)/(K_EV_AC+vesicleInput(1));
ATP_fraction_B0 = m_ATP_B0/(Km_ADNCYC_ATP+m_ATP_B0);
Pi_fraction_B0 = m_Pi_B0/(Km_PUNP1_Pi+m_Pi_B0);
uracil_fraction_B0 = m_uracil_B0/(Km_PYNP2r_uracil+m_uracil_B0);
dUTP_fraction_B0 = m_dUTP_B0/(Km_r1156_dUTP+m_dUTP_B0);
V_ADNCYC_target = J_B0_ref/(EV_activity_base*ATP_fraction_B0);
V_PDE1_target = J_B0_ref/(cAMP_ref/(Km_PDE1_cAMP+cAMP_ref));
V_NTD7_target = J_B0_ref/(AMP_ref/(Km_NTD7_AMP+AMP_ref));
V_PUNP1_target = J_B0_ref/((adn_ref/(Km_PUNP1_adn+adn_ref))*Pi_fraction_B0);
V_PYNP2r_target = J_B0_ref/((r1p_ref/(Km_PYNP2r_r1p+r1p_ref))*uracil_fraction_B0);
V_r1156_target = J_B0_ref/((uri_ref/(Km_r1156_uri+uri_ref))*dUTP_fraction_B0);
E_ADNCYC = V_ADNCYC_target/kcat_ADNCYC;
E_PDE1 = V_PDE1_target/kcat_PDE1; E_NTD7 = V_NTD7_target/kcat_NTD7;
E_PUNP1 = V_PUNP1_target/kcat_PUNP1; E_PYNP2r = V_PYNP2r_target/kcat_PYNP2r;
E_r1156 = V_r1156_target/kcat_r1156;
V_ADNCYC = kcat_ADNCYC*E_ADNCYC;
V_PDE1 = kcat_PDE1*E_PDE1; V_NTD7 = kcat_NTD7*E_NTD7;
V_PUNP1 = kcat_PUNP1*E_PUNP1; V_PYNP2r = kcat_PYNP2r*E_PYNP2r;
V_r1156 = kcat_r1156*E_r1156;
block0_parameters = table( ...
    {'ADNCYC (assumed EV control)';'PDE1';'NTD7';'PUNP1';'PYNP2r (reverse)';'r1156'}, ...
    [Km_ADNCYC_ATP;Km_PDE1_cAMP;Km_NTD7_AMP;Km_PUNP1_adn;Km_PYNP2r_r1p;Km_r1156_uri], ...
    [NaN;NaN;NaN;Km_PUNP1_Pi;Km_PYNP2r_uracil;Km_r1156_dUTP], ...
    [kcat_ADNCYC;kcat_PDE1;kcat_NTD7;kcat_PUNP1;kcat_PYNP2r;kcat_r1156], ...
    [E_ADNCYC;E_PDE1;E_NTD7;E_PUNP1;E_PYNP2r;E_r1156], ...
    [V_ADNCYC;V_PDE1;V_NTD7;V_PUNP1;V_PYNP2r;V_r1156], ...
    'VariableNames',{'Reaction','Km_first_mM','Km_second_mM', ...
                     'kcat_per_s','Enzyme_mM','Vmax_mM_per_s'});
disp(block0_parameters);

% Counter: buffered production activated ONLY by Y2.
E_Cp = 1e-3; SA_Cp = 10^0.5; MW_Cp = 75;
kcat_Cp = SA_Cp*MW_Cp/60; V_Cp = kcat_Cp*E_Cp;
E_Cr = 1e-3; SA_Cr = 10^1.5; MW_Cr = 75;
kcat_Cr = SA_Cr*MW_Cr/60; V_Cr = kcat_Cr*E_Cr;
Km_Cr_Cs = 10^0.5;
E_Cg = 1; SA_Cg = 10^1.5; MW_Cg = 275;
kcat_Cg = SA_Cg*MW_Cg/60; V_Cg = kcat_Cg*E_Cg;
Km_Cg_c = 10^0.5; Km_Cg_Cs = 10^0.5;
m_Cs_base = 0.01; k_Cs_supply = 1;

% Single affine neuron: two sources converge onto Z=m_Ni (paper Eq. (6)).
% r_w: C + buffered precursor -> C + Z, rate k_Nw*C.
% r_b: B + buffered precursor -> B + Z, rate k_Nb*B.
% r_I: Z -> products, rate k_Ni_deg*Z (single-pool removal).
% These catalytic source motifs preserve the upstream counter concentration.
% The linear laws are effective first-order kinetics; an unsaturated enzyme
% has k_eff = kcat*[E]/Km. Saturating the input changes the affine mapping.
w_N = 1;                           % dimensionless positive weight, illustrative
b_N = 1e-3;                        % mM, independent bias, illustrative
tau_Ni = 1;                        % s, common-pool relaxation time
NN_activation = 'identity';         % 'identity' (affine) or 'hill'
assert(isfinite(w_N) && w_N >= 0 && isfinite(b_N) && b_N >= 0, ...
    'This positive-concentration neuron requires w_N >= 0 and b_N >= 0.');
assert(isfinite(tau_Ni) && tau_Ni > 0, 'tau_Ni must be finite and positive.');
assert(strcmpi(NN_activation,'identity') || strcmpi(NN_activation,'hill'), ...
    'Choose NN_activation = identity or hill.');
k_Ni_deg = 1/tau_Ni;                % 1/s, removal rate k_r in Eq. (6)
k_Nw = k_Ni_deg*w_N;                % 1/s, k_Nw/k_Ni_deg = w_N
k_Nb = k_Ni_deg;                    % 1/s, bias coefficient is ONE
m_b_const = b_N;                    % mM, maintained bias input B
% w_N and b_N are not fitted distance-regression coefficients. m_c and
% m_out are concentrations, not time and distance; any decoding needs calibration.
% Optional output activation retains the previous Hill production/decay law.
V_No = 1; K_No = 0.1; n_No = 2;
k_out_deg = 1;
assert(V_No >= 0 && isfinite(V_No) && K_No > 0 && isfinite(K_No) && ...
    n_No >= 1 && isfinite(n_No) && k_out_deg > 0 && isfinite(k_out_deg), ...
    'Check optional Hill-output parameters.');

% Stage-2 integration controls: explicit Euler, with smaller local steps.
euler_safety_2 = 0.25;
max_substeps_2 = 20000;       % stop rather than hang in a very stiff regime
max_step_2 = 0.05;             % s; data sampling interval is unchanged
%%
%[text] ## Initial concentrations: equilibrium at the first input sample
root_options = optimset('Display','off','TolX',realmin);
% Equilibrate the adapter with the downstream consuming load.
J_B0_base = V_ADNCYC*EV_activity_base*ATP_fraction_B0;
assert(J_B0_base > 0 && J_B0_base < V_PDE1 && J_B0_base < V_NTD7 && ...
    J_B0_base < V_PUNP1*Pi_fraction_B0 && ...
    J_B0_base < V_PYNP2r*uracil_fraction_B0 && ...
    J_B0_base < V_r1156*dUTP_fraction_B0, ...
    'Block 0 has no finite equilibrium: source flux exceeds an adapter capacity.');
cAMP_base = Km_PDE1_cAMP*J_B0_base/(V_PDE1-J_B0_base);
AMP_base = Km_NTD7_AMP*J_B0_base/(V_NTD7-J_B0_base);
adn_base = Km_PUNP1_adn*J_B0_base/(V_PUNP1*Pi_fraction_B0-J_B0_base);
r1p_base = Km_PYNP2r_r1p*J_B0_base/(V_PYNP2r*uracil_fraction_B0-J_B0_base);
uri_base = Km_r1156_uri*J_B0_base/(V_r1156*dUTP_fraction_B0-J_B0_base);
% --- BRENDA-Km edit (3/4): the two UMP sources have DIFFERENT Km values ---
% Before (valid only for KM_F = KM_S):
%   input_fraction_base = J_B0_base/(V_S+V_F_D1);
%   u_ref = KM_S*input_fraction_base/(1-input_fraction_base);
%   j1_base = V_S*input_fraction_base;  jfast1_base = V_F_D1*input_fraction_base;
% Now the basal UMP is the root of  supply = r_S withdrawal + r_F withdrawal.
input_fraction_base = J_B0_base/(V_S+V_F);
assert(input_fraction_base > 0 && input_fraction_base < 1, ...
    'Block 0 supply exceeds the two UMP-source capacities.');
UMP_balance = @(u) V_S*u/(KM_S+u)+V_F*u/(KM_F+u)-J_B0_base;
u_upper = max(KM_S,KM_F)*input_fraction_base/(1-input_fraction_base);  % withdrawal >= supply here
u_ref = fzero(UMP_balance,[0 u_upper],root_options);
j1_base = V_S*u_ref/(KM_S+u_ref);
jfast1_base = V_F*u_ref/(KM_F+u_ref);
% --- end of edit (3/4) ---
assert(j1_base < V_Sr, 'No finite slow-memory equilibrium.');
M1_base = KM_Sr*j1_base/(V_Sr-j1_base);
% --- BRENDA-Km edit (6/7): equilibrium of F1 and Y1 with KM_I different from KM_Fr ---
% Before (valid only for KM_I = KM_Fr):
%   % Ignoring downstream load gives an upper bound for the Y1 root.
%   isolated_cancel_capacity = V_I_D1/V_Fr*(jfast1_base-j1_base);
%   assert(isolated_cancel_capacity > j1_base, ...
%       'No finite isolated D1 equilibrium; increase cancellation capacity.');
%   Y1_isolated = KM_I*j1_base/(isolated_cancel_capacity-j1_base);
%   balance_initial_Y1 = @(z) j1_base-load_initial_Y1(Y1_isolated*z) - ...
%       V_I_D1/V_Fr*(jfast1_base-j1_base+load_initial_Y1(Y1_isolated*z)) * ...
%       (Y1_isolated*z)/(KM_I+Y1_isolated*z);
%   Y1_base = Y1_isolated*fzero(balance_initial_Y1,[0 1],root_options);
% Now, for a trial Y1:  v_I = j1_base - load(Y1)             (Y1 balance)
%                       v_Fk = jfast1_base - v_I  ->  F1      (F1 balance)
% and Y1 is the root of  V_I*F1/(KM_I+F1)*Y1/(KM_I+Y1) = v_I.
load_initial_Y1 = @(y) consume_interstage*(V_F_D2+V_S1_D2)*y/(Km_D2_X+y);
vI_initial = @(y) j1_base-load_initial_Y1(y);
F1_initial = @(y) KM_Fr*(jfast1_base-vI_initial(y))/(V_Fr-(jfast1_base-vI_initial(y)));
balance_initial_Y1 = @(y) V_I*F1_initial(y)/(KM_I+F1_initial(y))*y/(KM_I+y)-vI_initial(y);
Y1_upper = 1e-3*KM_I;
% while true
%     % assert(V_Fr > jfast1_base-vI_initial(Y1_upper) && Y1_upper < 1e6*KM_I, ...
%     %     'No finite D1 equilibrium; increase cancellation capacity.');
%     if balance_initial_Y1(Y1_upper) >= 0, break; end
%     Y1_upper = 2*Y1_upper;
% end
Y1_base = 0;%fzero(balance_initial_Y1,[0 Y1_upper],root_options);
% --- end of edit (6/7) ---
j_I_D1_base = j1_base-load_initial_Y1(Y1_base);
j_R_D1_base = jfast1_base-j_I_D1_base;
% assert(j_I_D1_base > 0 && j_R_D1_base < V_Fr, ...
%     'No finite loaded fast-pool equilibrium.');
F1_delay_base = KM_Fr*j_R_D1_base/(V_Fr-j_R_D1_base);

j_F_D2_base = V_F_D2*Y1_base/(Km_D2_X+Y1_base);
j_S1_D2_base = V_S1_D2*Y1_base/(Km_D2_X+Y1_base);
% assert(V_R_D2 > j_F_D2_base && V_S2_D2 > j_S1_D2_base, ...
%     'Check stage-2 equilibrium capacities after editing enzyme amounts.');
% At equilibrium, v_I_D2=j_F_D2-v_R_D2 and v_S2_D2=j_S1_D2-j_F_D2+v_R_D2.
upper_Y2 = Km_R_D2_Y*j_F_D2_base/(V_R_D2-j_F_D2_base);
rate_R_initial = @(y) V_R_D2*y/(Km_R_D2_Y+y);
balance_initial_Y2 = @(z) j_F_D2_base - rate_R_initial(upper_Y2*z) - ...
    V_I_D2/V_S2_D2*(j_S1_D2_base-j_F_D2_base+rate_R_initial(upper_Y2*z)) * ...
    (upper_Y2*z)/(Km_I_D2_Y+upper_Y2*z);
Y2_base = upper_Y2*fzero(balance_initial_Y2,[0 1],root_options);
j_R_D2_base = rate_R_initial(Y2_base);
memory_fraction_base = (j_S1_D2_base-j_F_D2_base+j_R_D2_base)/V_S2_D2;
M2_base = Km_D2_M*memory_fraction_base/(1-memory_fraction_base);

gate_base = Y2_base/(K_gate2+Y2_base);
j_Cp_base = V_Cp*gate_base;
assert(j_Cp_base < k_Cs_supply*m_Cs_base, 'Counter source exceeds reset supply.');
balance_initial_Cs = @(z) k_Cs_supply*(m_Cs_base-m_Cs_base*z) - ...
    V_Cr*gate_base*(m_Cs_base*z)/(Km_Cr_Cs+m_Cs_base*z)-j_Cp_base;
Cs_base = m_Cs_base*fzero(balance_initial_Cs,[0 1],root_options);
reset_capacity_base = V_Cg*Cs_base/(Km_Cg_Cs+Cs_base);
assert(j_Cp_base < reset_capacity_base, 'Counter source exceeds reset capacity.');
counter_base = Km_Cg_c*j_Cp_base/(reset_capacity_base-j_Cp_base);
Na_base = k_Nw*counter_base/k_Ni_deg; % weighted contribution to Z
Nb_base = k_Nb*m_b_const/k_Ni_deg;   % independent bias contribution to Z
Ni_base = Na_base+Nb_base;          % common physical pool: w_N*C+b_N
if strcmpi(NN_activation,'hill')
    out_base = V_No*hill_activation(Ni_base,K_No,n_No)/k_out_deg;
else
    out_base = Ni_base;             % identity: output is the same physical pool
end
Y1_ref = Y1_base; Y2_ref = Y2_base;
%%
%[text] ## Concentration arrays
m_cAMP = zeros(N,1); m_cAMP(1) = cAMP_base;
m_AMP = zeros(N,1); m_AMP(1) = AMP_base;
m_adn = zeros(N,1); m_adn(1) = adn_base;
m_r1p = zeros(N,1); m_r1p(1) = r1p_base;
m_uri = zeros(N,1); m_uri(1) = uri_base;
m_in = zeros(N,1); m_in(1) = u_ref; % dynamic UMP, supplied only by r1156
v_ADNCYC = zeros(N,1); v_PDE1 = zeros(N,1); v_NTD7 = zeros(N,1);
v_PUNP1 = zeros(N,1); v_PYNP2r = zeros(N,1); v_r1156 = zeros(N,1);
m_S = zeros(N,1); m_S(1) = M1_base;
m_F = zeros(N,1); m_F(1) = F1_delay_base;
m_d1 = zeros(N,1); m_d1(1) = Y1_base;
m_M2 = zeros(N,1); m_M2(1) = M2_base;
m_Y2 = zeros(N,1); m_Y2(1) = Y2_base;
m_d1 = zeros(N,1); m_d1(1) = Y1_base;
m_d2 = zeros(N,1); m_d2(1) = Y2_base;
m_c = zeros(N,1); m_c(1) = counter_base;
m_Cs = zeros(N,1); m_Cs(1) = Cs_base;
m_Ni = zeros(N,1); m_Ni(1) = Ni_base;
m_Na = zeros(N,1); m_Na(1) = Na_base; % diagnostic weighted share of m_Ni
m_Nb = zeros(N,1); m_Nb(1) = Nb_base; % diagnostic bias share of m_Ni
m_out = zeros(N,1); m_out(1) = out_base;
v_S = zeros(N,1); v_Sk = zeros(N,1); v_F = zeros(N,1);
v_Fk = zeros(N,1); v_I_D1 = zeros(N,1);
v_F_D2_arr = zeros(N,1); v_S1_D2_arr = zeros(N,1); v_S2_D2_arr = zeros(N,1);
v_I_D2_arr = zeros(N,1); v_R_D2_arr = zeros(N,1);
v_Cp_arr = zeros(N,1); v_Cr_arr = zeros(N,1); v_Cg_arr = zeros(N,1);
enzyme_active = zeros(N,1);
%%
%[text] ## Linear reference parameters: plots only
d_S1_D1 = V_S*KM_S/(KM_S+u_ref)^2;
d_S2_D1 = V_Sr*KM_Sr/(KM_Sr+M1_base)^2;
d_F_D1_X = V_F*KM_F/(KM_F+u_ref)^2;
d_R_D1_F = V_Fr*KM_Fr/(KM_Fr+F1_delay_base)^2;
d_I_D1_F = V_I*Y1_base/(KM_I+Y1_base) * ...
    KM_I/(KM_I+F1_delay_base)^2;
d_I_D1_Y = V_I*F1_delay_base/(KM_I+F1_delay_base) * ...
    KM_I/(KM_I+Y1_base)^2;

d_F_D2 = V_F_D2*Km_D2_X/(Km_D2_X+Y1_base)^2;
d_S1_D2 = V_S1_D2*Km_D2_X/(Km_D2_X+Y1_base)^2;
d_S2_D2_M = V_S2_D2*Km_D2_M/(Km_D2_M+M2_base)^2;
d_I_D2_M = V_I_D2*Y2_base/(Km_I_D2_Y+Y2_base)*Km_D2_M/(Km_D2_M+M2_base)^2;
d_I_D2_Y = V_I_D2*M2_base/(Km_D2_M+M2_base)*Km_I_D2_Y/(Km_I_D2_Y+Y2_base)^2;
d_R_D2_Y = V_R_D2*Km_R_D2_Y/(Km_R_D2_Y+Y2_base)^2;
reference_det_2 = (d_S2_D2_M+d_I_D2_M)*(d_I_D2_Y+d_R_D2_Y)-d_I_D2_M*d_I_D2_Y;
positive_derivative_gain_2 = d_F_D2/reference_det_2;
% Loaded local response times, including the shared cancellation coupling.
% These are relaxation times, not fixed delays or prescribed peak times.
reference_trace_1 = d_R_D1_F+d_I_D1_F+d_I_D1_Y + ...
    consume_interstage*(d_F_D2+d_S1_D2);
reference_det_1 = (d_R_D1_F+d_I_D1_F) * ...
    (d_I_D1_Y+consume_interstage*(d_F_D2+d_S1_D2))-d_I_D1_F*d_I_D1_Y;
reference_root_1 = sqrt(reference_trace_1^2-4*reference_det_1);
tau_loaded_D1_memory = 1/d_S2_D1;
tau_loaded_D1_slow = (reference_trace_1+reference_root_1)/(2*reference_det_1);
tau_loaded_D1_fast = 2/(reference_trace_1+reference_root_1);
reference_trace_2 = d_S2_D2_M+d_I_D2_M+d_I_D2_Y+d_R_D2_Y;
reference_root_2 = sqrt(reference_trace_2^2-4*reference_det_2);
tau_loaded_D2_slow = (reference_trace_2+reference_root_2)/(2*reference_det_2);
tau_loaded_D2_fast = 2/(reference_trace_2+reference_root_2);
m_M1_ref = zeros(N,1); m_Y1_ref = zeros(N,1);
m_F1_delay_ref = zeros(N,1);
m_M2_ref = zeros(N,1); m_Y2_ref = zeros(N,1);
%%
%[text] ## Discrete-time simulation
for n = 1:N-1
    dt0 = t(n+1)-t(n);

    % ============================================================
    % BLOCK 0: ADNCYC anchor, then the five SELECTED adapter reactions
    % ============================================================
    % EV activates ATP -> cAMP + PPi; EV is not consumed by this law.
    v_ADNCYC(n) = V_ADNCYC*vesicleInput(n)/(K_EV_AC+vesicleInput(n)) * ...
        m_ATP_B0/(Km_ADNCYC_ATP+m_ATP_B0);
    % PDE1: cAMP + H2O -> AMP + H
    v_PDE1(n) = V_PDE1*m_cAMP(n)/(Km_PDE1_cAMP+m_cAMP(n));
    % NTD7: AMP + H2O -> adenosine + Pi
    v_NTD7(n) = V_NTD7*m_AMP(n)/(Km_NTD7_AMP+m_AMP(n));
    % PUNP1: adenosine + Pi -> adenine + ribose-1-phosphate
    v_PUNP1(n) = V_PUNP1*m_adn(n)/(Km_PUNP1_adn+m_adn(n)) * ...
        m_Pi_B0/(Km_PUNP1_Pi+m_Pi_B0);
    % PYNP2r REVERSE: ribose-1-phosphate + uracil -> uridine + Pi
    v_PYNP2r(n) = V_PYNP2r*m_r1p(n)/(Km_PYNP2r_r1p+m_r1p(n)) * ...
        m_uracil_B0/(Km_PYNP2r_uracil+m_uracil_B0);
    % r1156: uridine + dUTP -> UMP + dUDP + H
    v_r1156(n) = V_r1156*m_uri(n)/(Km_r1156_uri+m_uri(n)) * ...
        m_dUTP_B0/(Km_r1156_dUTP+m_dUTP_B0);

    % ============================================================
    % FIRST DIFFERENTIATOR: TWO LAG BRANCHES + SHARED CANCELLATION
    % ============================================================
    % Slow producer/consumer pair: UMP -> M1 -> Y1.
    v_S(n) = V_S*m_in(n)/(KM_S+m_in(n));
    v_Sk(n) = V_Sr*m_S(n)/(KM_Sr+m_S(n));
    % Fast producer/consumer pair: UMP -> F1_delay -> products.
    v_F(n) = V_F*m_in(n)/(KM_F+m_in(n));
    v_Fk(n) = V_Fr*m_F(n)/(KM_Fr+m_F(n));
    % r_I_D1 consumes BOTH the fast pool and m_d1, at the SAME flux.
    v_I_D1(n) = V_I*m_F(n)/(KM_I+m_F(n)) * ...
        m_d1(n)/(KM_I+m_d1(n));

    % The two second-stage sources withdraw Y1 in consuming mode.
    v_F_D2 = V_F_D2*m_d1(n)/(Km_D2_X+m_d1(n));
    v_S1_D2 = V_S1_D2*m_d1(n)/(Km_D2_X+m_d1(n));
    m_S(n+1) = m_S(n)+dt0*(v_S(n)-v_Sk(n));
    m_F(n+1) = m_F(n)+dt0*(v_F(n)-v_Fk(n)-v_I_D1(n));
    m_d1(n+1) = m_d1(n)+dt0*(v_Sk(n)-v_I_D1(n));
    % m_d1(n+1) = m_d1(n)+dt0*(v_Sk(n)-v_I_D1(n) - ...
    %     consume_interstage*(v_F_D2+v_S1_D2));


    % All Block 0 balances use rates at n, as ordinary explicit Euler.
    m_cAMP(n+1) = m_cAMP(n)+dt0*(v_ADNCYC(n)-v_PDE1(n));
    m_AMP(n+1) = m_AMP(n)+dt0*(v_PDE1(n)-v_NTD7(n));
    m_adn(n+1) = m_adn(n)+dt0*(v_NTD7(n)-v_PUNP1(n));
    m_r1p(n+1) = m_r1p(n)+dt0*(v_PUNP1(n)-v_PYNP2r(n));
    m_uri(n+1) = m_uri(n)+dt0*(v_PYNP2r(n)-v_r1156(n));
    % Actual differentiator input: UMP production minus BOTH UMP withdrawals.
    m_in(n+1) = m_in(n)+dt0*(v_r1156(n)-v_S(n)-v_F(n));


    % ============================================================
    % SECOND DIFFERENTIATOR: POSITIVE polarity and output cleanup
    % ============================================================
    % Fast source r_F_D2 produces Y2 directly. r_S1_D2 produces its inhibitor M2.
    % Hold Y1(n) fixed over this input interval, as in explicit Euler.
    M2_current = m_M2(n);
    Y2_current = m_Y2(n);
    time_left_2 = dt0;
    substeps_2 = 0;
    v_F_D2_arr(n) = v_F_D2; v_S1_D2_arr(n) = v_S1_D2;
    v_S2_D2_arr(n) = V_S2_D2*m_M2(n)/(Km_D2_M+m_M2(n));
    v_I_D2_arr(n) = V_I_D2*m_M2(n)/(Km_D2_M+m_M2(n))*m_Y2(n)/(Km_I_D2_Y+m_Y2(n));
    v_R_D2_arr(n) = V_R_D2*m_Y2(n)/(Km_R_D2_Y+m_Y2(n));
    % 
    % % Smaller Euler steps are needed when consuming inhibition is fast.
    % % These are ordinary explicit updates of the SAME two species.
    % while time_left_2 > 1e-12
    %     substeps_2 = substeps_2+1;
    %     % assert(substeps_2 <= max_substeps_2, ...
    %     %     ['Stage 2 is too stiff for this explicit Euler configuration at n=%d. ' ...
    %     %      'Its assumed Km_I_D2_Y can force extremely small steps. ' ...
    %     %      'Recalibrate the stage-2 kinetics or use stiff integration for this regime.'],n);
    %     % r_S2_D2: inhibitor M2 -> products
    %     v_S2_D2 = V_S2_D2*M2_current/(Km_D2_M+M2_current);
    %     % r_I_D2: M2 + Y2 -> products; SAME rate removes both species
    %     v_I_D2 = V_I_D2*M2_current/(Km_D2_M+M2_current) * ...
    %         Y2_current/(Km_I_D2_Y+Y2_current);
    %     % r_R_D2: extra Y2 removal for fast pulse cleanup
    %     v_R_D2 = V_R_D2*Y2_current/(Km_R_D2_Y+Y2_current);
    % 
    %     net_M2 = v_S1_D2-v_S2_D2-v_I_D2;
    %     net_Y2 = v_F_D2-v_I_D2-v_R_D2;
    % 
    %     % Limit the local Euler step using the current relaxation rates.
    %     relaxation_M2 = (V_S2_D2+V_I_D2*Y2_current/(Km_I_D2_Y+Y2_current)) * ...
    %         Km_D2_M/(Km_D2_M+M2_current)^2;
    %     relaxation_Y2 = V_I_D2*M2_current/(Km_D2_M+M2_current) * ...
    %         Km_I_D2_Y/(Km_I_D2_Y+Y2_current)^2 + ...
    %         V_R_D2*Km_R_D2_Y/(Km_R_D2_Y+Y2_current)^2;
    %     dt2 = min([time_left_2,max_step_2, ...
    %         euler_safety_2/max(relaxation_M2,relaxation_Y2)]);
    %     % A negative net rate cannot consume more than the available pool.
    %     if net_M2 < 0
    %         dt2 = min(dt2,0.8*M2_current/(-net_M2));
    %     end
    %     if net_Y2 < 0
    %         dt2 = min(dt2,0.8*Y2_current/(-net_Y2));
    %     end
    %     assert(isfinite(dt2) && dt2 > 0, 'Invalid stage-2 Euler substep.');
    % 
    %     % Direct integration: no implicit solver, subtraction gate, or clip.
    %     M2_current = M2_current+dt2*net_M2;
    %     Y2_current = Y2_current+dt2*net_Y2;
    %     time_left_2 = time_left_2-dt2;
    % end
    m_M2(n+1) = M2_current;
    m_Y2(n+1) = Y2_current;
    m_d2(n+1) = m_Y2(n+1);   % physical positive-pulse counter input

    % ============================================================
    % COUNTER: ACTIVATED ONLY BY Y2, ordinary n=1 kinetics
    % ============================================================
    active_fraction = m_d2(n)/(K_gate2+m_d2(n));
    enzyme_active(n) = active_fraction;
    v_Cp = V_Cp*active_fraction;
    v_Cr = V_Cr*active_fraction*m_Cs(n)/(Km_Cr_Cs+m_Cs(n));
    v_Cg = V_Cg*m_c(n)/(Km_Cg_c+m_c(n))*m_Cs(n)/(Km_Cg_Cs+m_Cs(n));
    v_Cs_supply = k_Cs_supply*(m_Cs_base-m_Cs(n));
    m_c(n+1) = m_c(n)+dt0*(v_Cp-v_Cg);
    m_Cs(n+1) = m_Cs(n)+dt0*(v_Cs_supply-v_Cr-v_Cg);
    v_Cp_arr(n) = v_Cp; v_Cr_arr(n) = v_Cr; v_Cg_arr(n) = v_Cg;

    % ============================================================
    % SINGLE NEURON: weighted sum with an independent constant bias
    % ============================================================
    % Sources r_w and r_b converge onto ONE physical pool m_Ni.
    % Euler integrates tau_Ni*dZ/dt+Z = w_N*m_c+b_N.
    assert(dt0*k_Ni_deg <= 1, ...
        'NN Euler step exceeds tau_Ni. Reduce dt0 or increase tau_Ni.');
    v_Na = k_Nw*m_c(n);             % r_w: weighted counter-regulated source
    v_Nb = k_Nb*m_b_const;          % r_b: constant source independent of m_c
    v_Ni = v_Na+v_Nb;              % both reactions produce the same Z
    m_Ni(n+1) = m_Ni(n)+dt0*(v_Ni-k_Ni_deg*m_Ni(n));
    % Bookkeeping for the original plot names, NOT two extra metabolites.
    % Their sum equals m_Ni; neither diagnostic drives any reaction.
    m_Na(n+1) = m_Na(n)+dt0*(v_Na-k_Ni_deg*m_Na(n));
    m_Nb(n+1) = m_Nb(n)+dt0*(v_Nb-k_Ni_deg*m_Nb(n));
    m_Nsum = m_Ni(n);               % old scalar name now refers to physical Z
    if strcmpi(NN_activation,'hill')
        assert(dt0*k_out_deg <= 1, ...
            'Hill-output Euler step too large. Reduce dt0 or k_out_deg.');
        v_No = V_No*m_Nsum^n_No/(K_No^n_No+m_Nsum^n_No);
        m_out(n+1) = m_out(n)+dt0*(v_No-k_out_deg*m_out(n));
    else
        m_out(n+1) = m_Ni(n+1);     % identity readout: no added reaction or lag
    end

    % assert(all(isfinite([m_cAMP(n+1),m_AMP(n+1),m_adn(n+1),m_r1p(n+1),m_uri(n+1),m_in(n+1),m_S(n+1),m_F(n+1),m_d1(n+1),m_M2(n+1),m_Y2(n+1), ...
    %     m_c(n+1),m_Cs(n+1),m_Ni(n+1),m_Na(n+1),m_Nb(n+1),m_out(n+1)])) && ...
    %     all([m_cAMP(n+1),m_AMP(n+1),m_adn(n+1),m_r1p(n+1),m_uri(n+1),m_in(n+1),m_S(n+1),m_F(n+1),m_d1(n+1),m_M2(n+1),m_Y2(n+1), ...
    %     m_c(n+1),m_Cs(n+1),m_Ni(n+1),m_Na(n+1),m_Nb(n+1),m_out(n+1)] >= 0), ...
    %     'Explicit Euler step invalid at n=%d. Reduce input timestep dt0.',n);

    % ============================================================
    % LINEAR REFERENCES: plotting only, include the changed stage-2 law
    % ============================================================
    input_deviation = m_in(n)-u_ref;
    m_M1_ref(n+1) = m_M1_ref(n)+dt0*(d_S1_D1*input_deviation-d_S2_D1*m_M1_ref(n));
    m_F1_delay_ref(n+1) = m_F1_delay_ref(n)+dt0*(d_F_D1_X*input_deviation - ...
        (d_R_D1_F+d_I_D1_F)*m_F1_delay_ref(n)-d_I_D1_Y*m_Y1_ref(n));
    m_Y1_ref(n+1) = m_Y1_ref(n)+dt0*(d_S2_D1*m_M1_ref(n) - ...
        d_I_D1_F*m_F1_delay_ref(n) - ...
        (d_I_D1_Y+consume_interstage*(d_F_D2+d_S1_D2))*m_Y1_ref(n));

    % The linear reference can also be fast, so use the same visible
    % substep concept rather than an implicit matrix or pool solve.
    M2_ref_current = m_M2_ref(n); Y2_ref_current = m_Y2_ref(n);
    time_left_ref = dt0;
    while time_left_ref > 1e-12
        dt_ref = min([time_left_ref,max_step_2, ...
            euler_safety_2/max(d_S2_D2_M+d_I_D2_M,d_I_D2_Y+d_R_D2_Y)]);
        net_M2_ref = d_S1_D2*m_Y1_ref(n)-(d_S2_D2_M+d_I_D2_M)*M2_ref_current-d_I_D2_Y*Y2_ref_current;
        net_Y2_ref = d_F_D2*m_Y1_ref(n)-d_I_D2_M*M2_ref_current-(d_I_D2_Y+d_R_D2_Y)*Y2_ref_current;
        M2_ref_current = M2_ref_current+dt_ref*net_M2_ref;
        Y2_ref_current = Y2_ref_current+dt_ref*net_Y2_ref;
        time_left_ref = time_left_ref-dt_ref;
    end
    m_M2_ref(n+1) = M2_ref_current; m_Y2_ref(n+1) = Y2_ref_current;
end

% Rates at the final sample, for inspecting each physical reaction.
v_ADNCYC(N) = V_ADNCYC*vesicleInput(N)/(K_EV_AC+vesicleInput(N))*ATP_fraction_B0;
v_PDE1(N) = V_PDE1*m_cAMP(N)/(Km_PDE1_cAMP+m_cAMP(N));
v_NTD7(N) = V_NTD7*m_AMP(N)/(Km_NTD7_AMP+m_AMP(N));
v_PUNP1(N) = V_PUNP1*m_adn(N)/(Km_PUNP1_adn+m_adn(N))*Pi_fraction_B0;
v_PYNP2r(N) = V_PYNP2r*m_r1p(N)/(Km_PYNP2r_r1p+m_r1p(N))*uracil_fraction_B0;
v_r1156(N) = V_r1156*m_uri(N)/(Km_r1156_uri+m_uri(N))*dUTP_fraction_B0;
v_S(N) = V_S*m_in(N)/(KM_S+m_in(N));
v_Sk(N) = V_Sr*m_S(N)/(KM_Sr+m_S(N));
v_F(N) = V_F*m_in(N)/(KM_F+m_in(N));
v_Fk(N) = V_Fr*m_F(N)/(KM_Fr+m_F(N));
v_I_D1(N) = V_I*m_F(N)/(KM_I+m_F(N))*m_d1(N)/(KM_I+m_d1(N));
v_F_D2_arr(N) = V_F_D2*m_d1(N)/(Km_D2_X+m_d1(N));
v_S1_D2_arr(N) = V_S1_D2*m_d1(N)/(Km_D2_X+m_d1(N));
v_S2_D2_arr(N) = V_S2_D2*m_M2(N)/(Km_D2_M+m_M2(N));
v_I_D2_arr(N) = V_I_D2*m_M2(N)/(Km_D2_M+m_M2(N))*m_Y2(N)/(Km_I_D2_Y+m_Y2(N));
v_R_D2_arr(N) = V_R_D2*m_Y2(N)/(Km_R_D2_Y+m_Y2(N));
enzyme_active(N) = m_Y2(N)/(K_gate2+m_Y2(N));
v_Cp_arr(N) = V_Cp*enzyme_active(N);
v_Cr_arr(N) = V_Cr*enzyme_active(N)*m_Cs(N)/(Km_Cr_Cs+m_Cs(N));
v_Cg_arr(N) = V_Cg*m_c(N)/(Km_Cg_c+m_c(N))*m_Cs(N)/(Km_Cg_Cs+m_Cs(N));
%%
%[text] ## Diagnostics: reference curves never enter a physical reaction rate
m_UMP = m_in; % alias only: m_in is the UMP pool integrated in the loop
input_consumption = v_S+v_F;
block0_inventory = m_cAMP+m_AMP+m_adn+m_r1p+m_uri+m_in;
block0_balance_residual = diff(block0_inventory) - ...
    diff(t).*(v_ADNCYC(1:end-1)-input_consumption(1:end-1));
fprintf('Block 0: basal cAMP %.4g mM; basal UMP m_in %.4g mM.\n',cAMP_base,u_ref);
fprintf('Block 0 tracked-pool balance residual %.4g mM per step.\n', ...
    max(abs(block0_balance_residual)));
m_Cs_arr = m_Cs;
d1_deviation = m_d1-Y1_ref;
d2_deviation = m_Y2-Y2_ref;
enable_concentration_ref = m_Y2; % physical output, no baseline subtraction
negativeSlope = -gradient(m_in,t);
risingSlope1 = gradient(m_d1,t); % diagnostic only; never drives a reaction
positive_slope2_ref = positive_derivative_gain_2*max(risingSlope1,0);
first_filtered_ref = m_Y1_ref;
second_filtered_enable_ref = Y2_ref+m_Y2_ref;
gate2_ref = gradient(m_d1,t) > 0; % digital rising-time diagnostic only
time_count_ref = [0; cumsum(diff(t).*double(gate2_ref(1:end-1)))];
counter_deviation = m_c-counter_base;
m_S1 = m_S; m_F1 = m_d1; m_S2 = m_M2; m_F2 = m_Y2;
m_d_target = first_filtered_ref;
m_d2_encoded = m_Y2;
fprintf('Stage-2 basal Y2 %.4g mM; peak %.4g mM; peak/basal %.4g.\n', ...
    Y2_ref,max(m_Y2),max(m_Y2)/Y2_ref);
fprintf('Counter activation: basal %.4g, peak %.4g.\n',gate_base,max(enzyme_active));

% Timing calibration: compare each output with the slope it should encode.
% A slope-transition duration differs from the time of its maximum slope.
% For a constant trace there is no meaningful slope-peak lag.
peak_slope_time_1_s = NaN; peak_output_time_1_s = NaN;
peak_slope_time_2_s = NaN; peak_output_time_2_s = NaN;
peak_lag_1_s = NaN; peak_lag_2_s = NaN;
if max(m_in)-min(m_in) > 100*eps*(1+max(abs(m_in)))
    [peak_slope_1,peak_slope_index_1] = max(negativeSlope);
    [peak_output_1,peak_output_index_1] = max(d1_deviation);
    if peak_slope_1 > 0 && peak_output_1 > 0
        peak_slope_time_1_s = t(peak_slope_index_1);
        peak_output_time_1_s = t(peak_output_index_1);
        peak_lag_1_s = peak_output_time_1_s-peak_slope_time_1_s;
    end
end
if max(m_d1)-min(m_d1) > 100*eps*(1+max(abs(m_d1)))
    [peak_slope_2,peak_slope_index_2] = max(risingSlope1);
    [peak_output_2,peak_output_index_2] = max(d2_deviation);
    if peak_slope_2 > 0 && peak_output_2 > 0
        peak_slope_time_2_s = t(peak_slope_index_2);
        peak_output_time_2_s = t(peak_output_index_2);
        peak_lag_2_s = peak_output_time_2_s-peak_slope_time_2_s;
    end
end
fprintf('Timing profile: %s. Time base and tau values: seconds.\n',timing_profile);
fprintf('D1 local relaxation times: memory %.3g s; coupled %.3g and %.3g s.\n', ...
    tau_loaded_D1_memory,tau_loaded_D1_slow,tau_loaded_D1_fast);
fprintf('D2 local relaxation times: %.3g and %.3g s.\n',tau_loaded_D2_slow,tau_loaded_D2_fast);
fprintf('D1: input slope peak %.3g s (%.3g min); output peak %.3g s (%.3g min); lag %.3g s.\n', ...
    peak_slope_time_1_s,peak_slope_time_1_s/60, ...
    peak_output_time_1_s,peak_output_time_1_s/60,peak_lag_1_s);
fprintf('D2: rising Y1 slope peak %.3g s (%.3g min); output peak %.3g s (%.3g min); lag %.3g s.\n', ...
    peak_slope_time_2_s,peak_slope_time_2_s/60, ...
    peak_output_time_2_s,peak_output_time_2_s/60,peak_lag_2_s);

%%
%[text] ## Plotting: original 13 figure roles, plus two Block 0 figures
shortWindow = [0,min(50,t(end))];
stage2Window = [0,min(300,t(end))];  % includes the Y1 peak and recovery
%%
%[text] ## UMP at the three plotted distances: Block 0 + its consuming first-stage load
%[text] Use the SAME kinetic constants and enzyme concentrations for every distance. The nine named front-end pools include the new fast memory; later blocks need not be rerun for these UMP traces.
idx = unique(round(linspace(1,size(vesicles,1),min(3,size(vesicles,1)))));
m_cAMP_dist = zeros(N,numel(idx)); m_AMP_dist = zeros(N,numel(idx));
m_adn_dist = zeros(N,numel(idx)); m_r1p_dist = zeros(N,numel(idx));
m_uri_dist = zeros(N,numel(idx)); m_in_dist = zeros(N,numel(idx));
m_M1_dist = zeros(N,numel(idx)); m_Y1_dist = zeros(N,numel(idx));
m_F1_delay_dist = zeros(N,numel(idx));

for j = 1:numel(idx)
    i_dist = idx(j);
    if i_dist == row
        % Reuse the exact selected-distance result from the main simulation.
        m_cAMP_dist(:,j) = m_cAMP; m_AMP_dist(:,j) = m_AMP;
        m_adn_dist(:,j) = m_adn; m_r1p_dist(:,j) = m_r1p;
        m_uri_dist(:,j) = m_uri; m_in_dist(:,j) = m_in;
        m_M1_dist(:,j) = m_S; m_Y1_dist(:,j) = m_d1;
        m_F1_delay_dist(:,j) = m_F;
        continue;
    end
    assert(all(isfinite(vesicles(i_dist,:))) && all(vesicles(i_dist,:) >= 0), ...
        'Each plotted distance requires a finite nonnegative EV trace.');

    % Equilibrium at this distance's first EV value, using common parameters.
    J_dist_base = V_ADNCYC*vesicles(i_dist,1)/(K_EV_AC+vesicles(i_dist,1))*ATP_fraction_B0;
    assert(J_dist_base < V_PDE1 && J_dist_base < V_NTD7 && ...
        J_dist_base < V_PUNP1*Pi_fraction_B0 && ...
        J_dist_base < V_PYNP2r*uracil_fraction_B0 && ...
        J_dist_base < V_r1156*dUTP_fraction_B0, ...
        'Adapter capacity exceeded at distance row %d.',i_dist);
    if J_dist_base > 0
        m_cAMP_dist(1,j) = Km_PDE1_cAMP*J_dist_base/(V_PDE1-J_dist_base);
        m_AMP_dist(1,j) = Km_NTD7_AMP*J_dist_base/(V_NTD7-J_dist_base);
        m_adn_dist(1,j) = Km_PUNP1_adn*J_dist_base/(V_PUNP1*Pi_fraction_B0-J_dist_base);
        m_r1p_dist(1,j) = Km_PYNP2r_r1p*J_dist_base/(V_PYNP2r*uracil_fraction_B0-J_dist_base);
        m_uri_dist(1,j) = Km_r1156_uri*J_dist_base/(V_r1156*dUTP_fraction_B0-J_dist_base);
        % --- BRENDA-Km edit (4/4): same root as in edit (3/4), for this distance ---
        % Before: m_in_dist(1,j) = KM_S*input_fraction_dist/(1-input_fraction_dist);
        %         J_S1_dist_base = V_S*input_fraction_dist;  J_F_dist_base = V_F_D1*input_fraction_dist;
        input_fraction_dist = J_dist_base/(V_S+V_F);
        assert(input_fraction_dist < 1, 'UMP load capacity exceeded at row %d.',i_dist);
        UMP_balance_dist = @(u) V_S*u/(KM_S+u)+V_F*u/(KM_F+u)-J_dist_base;
        u_upper_dist = max(KM_S,KM_F)*input_fraction_dist/(1-input_fraction_dist);
        m_in_dist(1,j) = fzero(UMP_balance_dist,[0 u_upper_dist],root_options);
        J_S1_dist_base = V_S*m_in_dist(1,j)/(KM_S+m_in_dist(1,j));
        J_F_dist_base = V_F*m_in_dist(1,j)/(KM_F+m_in_dist(1,j));
        % --- end of edit (4/4) ---
        assert(J_S1_dist_base < V_Sr, 'Slow memory capacity exceeded at row %d.',i_dist);
        m_M1_dist(1,j) = KM_Sr*J_S1_dist_base/(V_Sr-J_S1_dist_base);
        % --- BRENDA-Km edit (7/7): same equilibrium as edit (6/7), for this distance ---
        % Before (valid only for KM_I = KM_Fr):
        %   cancel_capacity_dist = V_I_D1/V_Fr*(J_F_dist_base-J_S1_dist_base);
        %   assert(cancel_capacity_dist > J_S1_dist_base, 'D1 cancellation capacity exceeded at row %d.',i_dist);
        %   upper_Y1_dist = KM_I*J_S1_dist_base/(cancel_capacity_dist-J_S1_dist_base);
        %   balance_Y1_dist = @(z) J_S1_dist_base-load_Y1_dist(upper_Y1_dist*z) - ...
        %   V_I_D1/V_Fr*(J_F_dist_base-J_S1_dist_base+load_Y1_dist(upper_Y1_dist*z)) * ...
        %   (upper_Y1_dist*z)/(KM_I+upper_Y1_dist*z);
        %   m_Y1_dist(1,j) = upper_Y1_dist*fzero(balance_Y1_dist,[0 1],root_options);
        load_Y1_dist = @(y) consume_interstage*(V_F_D2+V_S1_D2)*y/(Km_D2_X+y);
        vI_dist0 = @(y) J_S1_dist_base-load_Y1_dist(y);
        F1_dist0 = @(y) KM_Fr*(J_F_dist_base-vI_dist0(y))/(V_Fr-(J_F_dist_base-vI_dist0(y)));
        balance_Y1_dist = @(y) V_I*F1_dist0(y)/(KM_I+F1_dist0(y))*y/(KM_I+y)-vI_dist0(y);
        upper_Y1_dist = 1e-3*KM_I;
        % while true
        %     assert(V_Fr > J_F_dist_base-vI_dist0(upper_Y1_dist) && upper_Y1_dist < 1e6*KM_I, ...
        %         'D1 cancellation capacity exceeded at row %d.',i_dist);
        %     if balance_Y1_dist(upper_Y1_dist) >= 0, break; end
        %     upper_Y1_dist = 2*upper_Y1_dist;
        % end
        m_Y1_dist(1,j) = 0;% fzero(balance_Y1_dist,[0 upper_Y1_dist],root_options);
        % --- end of edit (7/7) ---
        J_R_dist_base = J_F_dist_base-J_S1_dist_base+load_Y1_dist(m_Y1_dist(1,j));
        % assert(J_R_dist_base < V_Fr, 'Fast memory capacity exceeded at row %d.',i_dist);
        m_F1_delay_dist(1,j) = KM_Fr*J_R_dist_base/(V_Fr-J_R_dist_base);

    end
    % For a zero initial EV signal, the zero-initialized pools are equilibrium.
    for n_dist = 1:N-1
        dt_dist = t(n_dist+1)-t(n_dist);
        v_AC_dist = V_ADNCYC*vesicles(i_dist,n_dist)/(K_EV_AC+vesicles(i_dist,n_dist))*ATP_fraction_B0;
        v_PDE_dist = V_PDE1*m_cAMP_dist(n_dist,j)/(Km_PDE1_cAMP+m_cAMP_dist(n_dist,j));
        v_NTD_dist = V_NTD7*m_AMP_dist(n_dist,j)/(Km_NTD7_AMP+m_AMP_dist(n_dist,j));
        v_PUN_dist = V_PUNP1*m_adn_dist(n_dist,j)/(Km_PUNP1_adn+m_adn_dist(n_dist,j))*Pi_fraction_B0;
        v_PYN_dist = V_PYNP2r*m_r1p_dist(n_dist,j)/(Km_PYNP2r_r1p+m_r1p_dist(n_dist,j))*uracil_fraction_B0;
        v_1156_dist = V_r1156*m_uri_dist(n_dist,j)/(Km_r1156_uri+m_uri_dist(n_dist,j))*dUTP_fraction_B0;
        v_S1_dist = V_S*m_in_dist(n_dist,j)/(KM_S+m_in_dist(n_dist,j));
        v_S2_dist = V_Sr*m_M1_dist(n_dist,j)/(KM_Sr+m_M1_dist(n_dist,j));
        v_F_dist = V_F*m_in_dist(n_dist,j)/(KM_F+m_in_dist(n_dist,j));
        v_R_dist = V_Fr*m_F1_delay_dist(n_dist,j)/(KM_Fr+m_F1_delay_dist(n_dist,j));
        v_I_dist = V_I*m_F1_delay_dist(n_dist,j)/(KM_I+m_F1_delay_dist(n_dist,j)) * ...
            m_Y1_dist(n_dist,j)/(KM_I+m_Y1_dist(n_dist,j));
        v_D2_load_dist = consume_interstage*(V_F_D2+V_S1_D2) * ...
            m_Y1_dist(n_dist,j)/(Km_D2_X+m_Y1_dist(n_dist,j));

        m_cAMP_dist(n_dist+1,j) = m_cAMP_dist(n_dist,j)+dt_dist*(v_AC_dist-v_PDE_dist);
        m_AMP_dist(n_dist+1,j) = m_AMP_dist(n_dist,j)+dt_dist*(v_PDE_dist-v_NTD_dist);
        m_adn_dist(n_dist+1,j) = m_adn_dist(n_dist,j)+dt_dist*(v_NTD_dist-v_PUN_dist);
        m_r1p_dist(n_dist+1,j) = m_r1p_dist(n_dist,j)+dt_dist*(v_PUN_dist-v_PYN_dist);
        m_uri_dist(n_dist+1,j) = m_uri_dist(n_dist,j)+dt_dist*(v_PYN_dist-v_1156_dist);
        m_in_dist(n_dist+1,j) = m_in_dist(n_dist,j)+dt_dist*(v_1156_dist-v_S1_dist-v_F_dist);
        m_M1_dist(n_dist+1,j) = m_M1_dist(n_dist,j)+dt_dist*(v_S1_dist-v_S2_dist);
        m_F1_delay_dist(n_dist+1,j) = m_F1_delay_dist(n_dist,j)+dt_dist*(v_F_dist-v_R_dist-v_I_dist);
        m_Y1_dist(n_dist+1,j) = m_Y1_dist(n_dist,j)+dt_dist*(v_S2_dist-v_I_dist-v_D2_load_dist);
        % assert(all(isfinite([m_cAMP_dist(n_dist+1,j),m_AMP_dist(n_dist+1,j), ...
        %     m_adn_dist(n_dist+1,j),m_r1p_dist(n_dist+1,j),m_uri_dist(n_dist+1,j), ...
        %     m_in_dist(n_dist+1,j),m_M1_dist(n_dist+1,j),m_F1_delay_dist(n_dist+1,j),m_Y1_dist(n_dist+1,j)])) && ...
        %     all([m_cAMP_dist(n_dist+1,j),m_AMP_dist(n_dist+1,j), ...
        %     m_adn_dist(n_dist+1,j),m_r1p_dist(n_dist+1,j),m_uri_dist(n_dist+1,j), ...
        %     m_in_dist(n_dist+1,j),m_M1_dist(n_dist+1,j),m_F1_delay_dist(n_dist+1,j),m_Y1_dist(n_dist+1,j)] >= 0), ...
        %     'Invalid Euler step at distance row %d; reduce dt0.',i_dist);
    end
end
slope_UMP_dist = diff(m_in_dist,1,1)./diff(t); % [mM/s], same order as idx

%%

% 1. Simulated UMP traces at three distances, with their finite-difference slopes.
figure;
plot(t/60,m_in_dist,'LineWidth',2); grid on;
xlabel('Time [min]'); ylabel('UMP: m_{in} [mM]');
labels = arrayfun(@(i) sprintf('d = %.3g um',dist(i)*1e6),idx,'UniformOutput',false);
legend(labels,'Location','southwest','NumColumns',numel(idx));
set(gca,'FontSize',fontsize);
positionMain = get(gca,'Position');
inset = axes('Position',[positionMain(1)+0.43 positionMain(2)+0.42 0.30 0.26]);
total = min(500,N-1);
plot(t(1:total),slope_UMP_dist(1:total,:),'LineWidth',1.5);
xlabel('Time [s]'); ylabel('UMP slope [mM/s]'); grid on;
set(inset,'FontSize',max(8,fontsize-2));

% 2. Actual input concentration: UMP at the output of Block 0.
figure;
plot(t/60,m_in,'LineWidth',2); grid on;
xlabel('Time [min]'); ylabel('UMP: m_{in} [mM]');
legend(sprintf('d = %.3g um',dist(row)*1e6),'Location','southwest');
set(gca,'FontSize',fontsize+4);

% 3. First-stage physical output, short time window, with its baseline.
figure;
plot(t/60,m_d1,'LineWidth',2); hold on; grid on;
xlabel('Time [min]');
ylabel('$[m_{d1}]$ [mM]','Interpreter','latex');
xlim([0,min(150,t(end)/60)]);
set(gca,'FontSize',fontsize+4);


% 4. Stage-2 physical output: same MINUTE axis/window as the first stage.
figure;
plot(t/60,m_Y2,'LineWidth',2); grid on;
ylabel('$[m_{d2}]$ [mM]','Interpreter','latex');
xlabel('Time [min]');
xlim([0,min(5,t(end)/60)]);
set(gca,'FontSize',fontsize+4);

% 5. Stage-2 physical output: same MINUTE axis/window as the first stage.
figure;
plot(t/60,m_c,'LineWidth',2); grid on;
ylabel('$[m_{c}]$ mM','Interpreter','latex');
xlabel('Time [min]');
xlim([0,min(5,t(end)/60)]);
set(gca,'FontSize',fontsize+4);

% yyaxis right; plot(t,time_count_ref,'--','LineWidth',2);
% ylabel('Digital Y_1 rising-time reference [s]');
% xlabel('Time [s]'); xlim(stage2Window); grid on;

% % 7. Physical Y2 and positive diagnostic: distinguish baseline and sign.
% figure;
% subplot(2,1,1);
% yyaxis left; plot(t,m_d1,'LineWidth',1.5); ylabel('Physical Y_1 [mM]');
% yyaxis right; plot(t,m_Y2,'LineWidth',1.5); hold on;
% yline(Y2_ref,'--'); ylabel('Physical Y_2 [mM]');
% xlabel('Time [s]'); xlim(stage2Window); grid on;
% legend('Y_1','Y_2','Y_{2,0}','Location','best');
% title('Rising Y_1 produces a positive Y_2 pulse');
% subplot(2,1,2);
% plot(t,enable_concentration_ref,'b', ...
%     t,second_filtered_enable_ref,'r--', ...
%     t,positive_slope2_ref,'k:','LineWidth',1.5);
% xlabel('Time [s]'); ylabel('Y_2 concentration [mM]');
% legend('Physical Y_2 output', ...
%     'Linearized-MM reference','Scaled positive digital slope reference', ...
%     'Location','best');
% xlim(stage2Window); grid on;
% title('Physical pulse and independent reference curves');
% 
% % 8. First-stage deviation against digital input slope.
% figure;
% yyaxis left; plot(t,d1_deviation,'LineWidth',2);
% ylabel('Y_1-Y_{1,0}, diagnostic [mM]');
% yyaxis right; plot(t,negativeSlope,'--','LineWidth',2);
% ylabel('-dm_{in}/dt [mM/s]'); xlabel('Time [s]');
% legend('Physical output deviation','Digital slope reference');
% xlim(shortWindow); grid on;
% 
% % 9. First-stage output deviation and its independent filtered reference.
% figure;
% plot(t,d1_deviation,t,first_filtered_ref,'--','LineWidth',2);
% legend('Physical Y_1-Y_{1,0}','Linearized coupled-MM reference');
% xlabel('Time [s]'); ylabel('Output deviation [mM]'); xlim(shortWindow); grid on;
% 
% % 10. Memory/output and received input.
% figure;
% yyaxis left; plot(t,m_S,t,m_d1,'LineWidth',2);
% ylabel('M_1, Y_1 [mM]');
% yyaxis right; plot(t,m_in,'LineWidth',2); ylabel('Input [mM]');
% xlabel('Time [s]'); legend('M_1','Y_1','Input','Location','best'); grid on;
% set(gca,'FontSize',fontsize);
% 
% % 11. Received input and first differentiator physical output.
% figure;
% yyaxis left; plot(t,m_in,'LineWidth',2); ylabel('Input [mM]');
% yyaxis right; plot(t,m_d1,'LineWidth',2); ylabel('Y_1 [mM]');
% xlabel('Time [s]'); legend('Input','First-stage output','Location','best'); grid on;
% title('First differentiator: output around a positive baseline');
% set(gca,'FontSize',fontsize);
% 
% % 12. Chemical counter output; show the basal level explicitly.
% figure;
% yyaxis left; plot(t,m_in,'LineWidth',2); ylabel('Input [mM]');
% yyaxis right; plot(t,m_c,'LineWidth',2); hold on;
% yline(counter_base,'--','Basal counter'); ylabel('Counter C [mM]');
% xlabel('Time [s]'); legend('Input','Chemical C','Basal C','Location','best'); grid on;
% title('Baseline-regulated analog counter'); set(gca,'FontSize',fontsize);
% 
% % 13. Affine NN response; original signal names and commented status retained.
% figure;
% plot(t,[m_c,m_Ni,m_Na,m_Nb,m_out],'LineWidth',2); grid on;
% xlabel('Time [s]'); ylabel('Concentration [mM]');
% legend('C: NN input','Z: common affine pool','Weighted share (diagnostic)', ...
%     'Bias share (diagnostic)','NN output','Location','best');
% title('NN block: w C + b with common-pool dynamics'); set(gca,'FontSize',fontsize);
% 
% 
% % 14. EV driver, generated cAMP, and actual UMP input; old surrogate is diagnostic only.
% figure;
% subplot(2,1,1);
% yyaxis left; plot(t/60,vesicleInput,'LineWidth',1.5); ylabel('EV signal [dataset units]');
% yyaxis right; plot(t/60,[m_cAMP,m_in],'LineWidth',1.5); ylabel('Concentration [mM]');
% xlabel('Time [min]'); grid on;
% legend('EV driver','Generated cAMP','UMP: actual m_{in}','Location','best');
% title('Block 0: assumed EV transduction and selected intracellular adapter');
% subplot(2,1,2);
% plot(t/60,m_in,t/60,legacy_input_ref,'--','LineWidth',1.5); grid on;
% xlabel('Time [min]'); ylabel('Concentration [mM]');
% legend('Dynamic UMP: actual m_{in}','Old direct EV scaling: comparison only','Location','best');
% title('The direct vesicle-to-input assignment no longer drives the reactions');
% 
% % 15. Adapter intermediates and UMP production/consumption balance.
% figure;
% subplot(2,1,1);
% plot(t/60,[m_cAMP,m_AMP,m_adn,m_r1p,m_uri,m_in],'LineWidth',1.5); grid on;
% xlabel('Time [min]'); ylabel('Concentration [mM]');
% legend('cAMP','AMP','Adenosine','Ribose-1-P','Uridine','UMP: m_{in}', ...
%     'Location','best','NumColumns',2);
% title('Block 0 pools (chosen equilibrium initialization)');
% subplot(2,1,2);
% plot(t/60,v_r1156,t/60,input_consumption, ...
%      t/60,v_r1156-input_consumption,'--','LineWidth',1.5); grid on;
% xlabel('Time [min]'); ylabel('Rate [mM/s]');
% legend('UMP production: r1156','UMP use: r_{S1,D1}+r_{F,D1}', ...
%     'Net dm_{in}/dt','Location','best');
% title('A declining UMP pool occurs when its consumption exceeds production');

% 16. Timing checks; the previously commented figures remain commented.
if plot_timing_checks
    figure;
    timingWindow = [0,min(5,t(end)/60)];
    subplot(2,1,1);
    yyaxis left;
    plot(t/60,d1_deviation,'LineWidth',1.5);
    ylabel('m_{d1}-Y_{1,0} [mM]');
    yyaxis right;
    plot(t/60,negativeSlope,'--','LineWidth',1.5);
    ylabel('-dm_{in}/dt [mM/s]');
    xlabel('Time [min]'); xlim(timingWindow); grid on;
    legend('First-stage output deviation','Input falling-slope reference','Location','best');
    title(sprintf('First-stage peak lag: %.3g s',peak_lag_1_s));

    subplot(2,1,2);
    yyaxis left;
    plot(t/60,m_d2,'LineWidth',1.5); hold on;
    yline(Y2_ref,':'); ylabel('m_{d2} [mM]');
    yyaxis right;
    plot(t/60,max(risingSlope1,0),'--','LineWidth',1.5);
    ylabel('Positive dm_{d1}/dt [mM/s]');
    xlabel('Time [min]'); xlim(timingWindow); grid on;
    legend('Second-stage output','Y_{2,0}','Y1 rising-slope reference','Location','best');
    title(sprintf('Second-stage peak lag: %.3g s',peak_lag_2_s));
end

%%
%[text] ## Helpers
function h = hill_activation(x,K,n)
    % Stable Hill fraction x^n/(K^n+x^n); zero input gives zero activation.
    h = 1./(1+exp(-n*log(x./K)));
end

function f = mm_fraction(x,Km)
    f = x./(Km+x);
end

function addingPathParentFolderByName(targetName)
    % Start from the current directory
    currFolder = pwd;
    found = false;
    
    % Continue searching until you reach the root folder
    while true
        % Get the parent folder
        [parentFolder, currentName] = fileparts(currFolder);

        % Check if the current folder's name is the target
        if strcmpi(currentName, targetName)
            found = true;
            break;
        end
        
        % If we've reached the root or no change, exit the loop
        if isempty(parentFolder) || strcmp(currFolder, parentFolder)
            break;
        end
        
        % Move one level up
        currFolder = parentFolder;
    end

    if found
        addpath(genpath(currFolder));
        fprintf('Adding matlab path to: %s\n', currFolder);
    else
        error('Folder named "%s" not found in any parent directory.', targetName);
    end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:6f7ea38d]
%   data: {"dataType":"text","outputData":{"text":"Adding matlab path to: D:\\code\\bioelectrodynamics\\code\n","truncated":false}}
%---
%[output:1f520a66]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Adapter 6: removal reaction \"r_S1\" not found."}}
%---
%[output:672d5d47]
%   data: {"dataType":"warning","outputData":{"text":"Warning: r_S: removal reaction \"r_S2\" not found."}}
%---
%[output:808e23bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: r_F: removal reaction \"r_R\" not found."}}
%---
%[output:770b142c]
%   data: {"dataType":"text","outputData":{"text":"                   <strong>Reaction<\/strong>         <strong>KM<\/strong>         <strong>kcat<\/strong>         <strong>E<\/strong>          <strong>k<\/strong>               <strong>RemovedBy<\/strong>               <strong>G<\/strong>            <strong>tau<\/strong>    \n                 <strong>_____________<\/strong>    <strong>_______<\/strong>    <strong>_________<\/strong>    <strong>_____<\/strong>    <strong>__________<\/strong>    <strong>______________________<\/strong>    <strong>__________<\/strong>    <strong>__________<\/strong>\n\n    <strong>Adapter 1<\/strong>    {'Adapter 1'}    0.00285       1.4608    0.001       0.51256    {'Adapter 2'         }       0.13904       0.27126\n    <strong>Adapter 2<\/strong>    {'Adapter 2'}     0.0245       90.318    0.001        3.6865    {'Adapter 3'         }        4.4025        1.1942\n    <strong>Adapter 3<\/strong>    {'Adapter 3'}      0.083         69.5    0.001       0.83735    {'Adapter 4'         }        2.4069        2.8745\n    <strong>Adapter 4<\/strong>    {'Adapter 4'}       0.44       153.07    0.001       0.34789    {'Adapter 5'         }        15.116        43.451\n    <strong>Adapter 5<\/strong>    {'Adapter 5'}      0.505       11.622    0.001      0.023014    {'Adapter 6, r_No'   }     0.0091461       0.39741\n    <strong>Adapter 6<\/strong>    {'Adapter 6'}      0.146          9.5    0.001      0.065068    {'r_S1, r_F'         }           NaN           NaN\n    <strong>r_S      <\/strong>    {'r_S'      }     0.0125       14.867    0.001        1.1893    {'r_S2'              }           NaN           NaN\n    <strong>r_Sr     <\/strong>    {'r_Sr'     }       0.06       126.03    0.001        2.1006    {'r_I*, r_S1_2, r_F2'}        1.9534       0.92995\n    <strong>r_F      <\/strong>    {'r_F'      }       0.36       11.622    0.001      0.032284    {'r_I, r_R'          }           NaN           NaN\n    <strong>r_I      <\/strong>    {'r_I'      }       0.33       126.03    0.001       0.38192    {0×0 char            }           NaN           NaN\n    <strong>r_Fr     <\/strong>    {'r_Fr'     }       0.02      0.26777    0.001      0.013388    {0×0 char            }           NaN           NaN\n    <strong>r_S1_2   <\/strong>    {'r_S1_2'   }       0.36       124.81    0.001        0.3467    {'r_S2_2, r_I2'      }    2.1769e-05    6.2789e-05\n    <strong>r_S2_2   <\/strong>    {'r_S2_2'   }      0.027      4.3e+05    0.001         15926    {'r_I2*, r_R2'       }         16356         1.027\n    <strong>r_F2     <\/strong>    {'r_F2'     }       0.36       124.81    0.001        0.3467    {'r_I2*, r_R2'       }       0.35607         1.027\n    <strong>r_I2     <\/strong>    {'r_I2'     }       0.36       124.81    0.001        0.3467    {0×0 char            }           NaN           NaN\n    <strong>r_R2     <\/strong>    {'r_R2'     }       0.42       263.33    0.001       0.62698    {0×0 char            }           NaN           NaN\n    <strong>r_Cp     <\/strong>    {'r_Cp'     }        1.3       3.2516    0.001     0.0025012    {'r_Cg'              }         5.348        2138.2\n    <strong>r_Cg     <\/strong>    {'r_Cg'     }        1.3        0.608    0.001    0.00046769    {0×0 char            }           NaN           NaN\n    <strong>r_Cr     <\/strong>    {'r_Cr'     }        1.8    0.0052112    0.001    2.8951e-06    {0×0 char            }           NaN           NaN\n    <strong>r_Ni     <\/strong>    {'r_Ni'     }      0.046      0.85091    0.001      0.018498    {'r_NiDecay'         }       0.14477        7.8262\n    <strong>r_NiDecay<\/strong>    {'r_NiDecay'}       0.01       1.2778    0.001       0.12778    {0×0 char            }           NaN           NaN\n    <strong>r_Na     <\/strong>    {'r_Na'     }    0.09465       8.0981    0.001      0.085559    {'r_NaDecay'         }        447.24        5227.3\n    <strong>r_NaDecay<\/strong>    {'r_NaDecay'}       11.5          2.2    0.001     0.0001913    {0×0 char            }           NaN           NaN\n    <strong>r_Nb     <\/strong>    {'r_Nb'     }    0.09465       8.0981    0.001      0.085559    {'r_NbDecay'         }        447.24        5227.3\n    <strong>r_NbDecay<\/strong>    {'r_NbDecay'}       11.5          2.2    0.001     0.0001913    {0×0 char            }           NaN           NaN\n    <strong>r_No     <\/strong>    {'r_No'     }     0.1595       390.97    0.001        2.4512    {'r_NoDecay'         }       0.95999       0.39163\n    <strong>r_NoDecay<\/strong>    {'r_NoDecay'}       0.55       1404.4    0.001        2.5534    {0×0 char            }           NaN           NaN\n\n","truncated":false}}
%---
%[output:3c07477a]
%   data: {"dataType":"text","outputData":{"text":"Algorithm 1: tau_S = 5.19 s, tau_F = 4.19 s, dt = 1 s, alpha = 0.807, WX*tau_S = 0.0865\n  [E]_S  =    0.1621 uM   k_S  = 0.1928 1\/s\n  [E]_Sr =   0.09177 uM   k_Sr = 0.1928 1\/s   (gain w_S = 1)\n  [E]_F  =     7.397 uM   k_F  = 0.2388 1\/s\n  [E]_Fr =     17.84 uM   k_Fr = 0.2388 1\/s   (gain w_F = 1, k_I = 0 1\/s)\n  feasible below 20 uM: 1;  linear regime: 1\n","truncated":false}}
%---
%[output:4bab2dfa]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"5.1875"}}
%---
%[output:2358b501]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"4.1875"}}
%---
%[output:0c8c3c4e]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.8072"}}
%---
%[output:2e69004f]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Unrecognized function or variable 'kcat_F_D2'."}}
%---
