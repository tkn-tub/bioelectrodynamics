clc;
% clearvars;
close all;

addingPathParentFolderByName('code'); %[output:4c4d511c]
Parameters                               % provides dt0 and fontsize
%%
%[text] 1\. Input data
load Dataset_cell2cell.mat
Dataset  = rmmissing(Dataset);
dist     = table2array(Dataset(:,1));
vesicles = table2array(Dataset(:,2:end));

row = 1;                                 % distance analysed in detail
N = size(vesicles,2);
t = (0:N-1)'*dt0;
assert(N >= 3 && all(isfinite(vesicles(row,:))) && all(vesicles(row,:) >= 0) && vesicles(row,1) > 0, ...
    'Use a finite nonnegative EV signal with a positive first sample.');
%%
%[text] Enzyme constants from the BRENDA workbook
%[text] Description: This block uses your credentials at "Brenda\\\_credentials.txt" for searching at \[Brenda website\](https://www.brenda-enzymes.info/functional\_parameter\_statistics.php) the enzyme parameters. Your credentials should be included in a txt file with two lines, the first line includes your email, and the second line your password. The enzymes are identified from the file "Block\\\_7\\\_Selected\\\_Motifs.txt", which is completed by the second step and includes the chemical reaction network of the algorithm. After the code reads the enzyme parameters it prepares an excel table with the listed values per enzyme.
% Uncomment below for reading the constants from https://www.brenda-enzymes.info/functional_parameter_statistics.php 

% %Loggin in Brenda for quering the enzyme parameters https://www.brenda-enzymes.info/functional_parameter_statistics.php 
% % Default values
% % Look for the credentials file on the MATLAB search path
% credentialsFile = which('BRENDA_credentials.txt');
% 
% if ~isempty(credentialsFile) % which returns empty if not found
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
%     password = passwordFromFile; 
% 
%     clear fid emailFromFile passwordFromFile
% else
%     error('BRENDA_credentials.txt not found on the MATLAB path.');
% end
% clear brendaKineticsTable
% brendaKineticsTable('Block_7_Selected_Motifs.txt', ...
%                    'kinetics_BRENDA.xlsx', email, password);

% Sheet "values" of kinetics_BRENDA.xlsx: reaction label in column B, and the
% columns K_M [mM], SA [umol/min/mg] and MW [mg/umol] found by their titles.
% kcat [1/s] = SA*MW/60.
raw = readcell('kinetics_BRENDA.xlsx','Sheet','values');
isText = @(x) ischar(x) || (isstring(x) && isscalar(x) && ~ismissing(x));
hdrRow = find(any(cellfun(@(x) isText(x) && startsWith(string(x),"K_M"),raw),2),1);
assert(~isempty(hdrRow),'Header row (a cell starting with "K_M") not found in the workbook.');
colOf = @(title) find(cellfun(@(x) isText(x) && startsWith(string(x),title),raw(hdrRow,:)),1);
colKM = colOf("K_M");  colSA = colOf("SA");  colMW = colOf("MW"); colKcat = colOf("k_cat");

%==================
% Block 0
%==================

% ---- Constants from the workbook ----
% Labels and chain exactly as the sheets "values" and "Reaction mapping" list them:
%   r_adapter_1  3',5'-cyclic AMP            -> AMP                  (K_M for cAMP)
%   r_adapter_2  AMP + argininosuccinate+PPi -> ATP                  (K_M for ATP)
%   r_adapter_3  ATP + aspartate + glutamine -> L-glutamate (m_d1)   (K_M for ATP)
%   r_adapter_4  L-glutamate + 3-P-hydroxypyruvate -> 3-P-L-serine   (K_M for glutamate)
%   r_adapter_5  3-P-L-serine                -> L-serine (m_in)      (K_M for 3-P-serine)
reactionNamesB0 = {'r_adapter_1','r_adapter_2','r_adapter_3','r_adapter_4'};
KM_B0 = nan(1,5);
kcat_B0 = nan(1,5);

for i = 1:5
    r = find(cellfun(@(x) isText(x) && strcmp(strtrim(char(x)),reactionNames{i}),raw(:,2)),1);
    assert(~isempty(r),'Reaction %s not found in column B of the workbook.',reactionNames{i});
    KM_B0(i)   = raw{r,colKM};
    %kcat(i) = raw{r,colSA}*raw{r,colMW}/60;
    kcat_B0(i) = raw{r,colKcat};
end
assert(all(isfinite([KM_B0 kcat_B0])) && all([KM_B0 kcat_B0] > 0),'K_M, SA and MW of the five reactions must be positive numbers.');


Km_a1_cAMP = KM_B0(1);  kcat_a1 = kcat_B0(1);   % cAMP        -> AMP
Km_a2_AMP  = KM_B0(2);  kcat_a2 = kcat_B0(2);   % AMP         -> ATP
Km_a3_ATP  = KM_B0(3);  kcat_a3 = kcat_B0(3);   % ATP         -> L-glutamate
Km_a4_glu  = KM_B0(4);  kcat_a4 = kcat_B0(4);   % L-glutamate -> 3-P-L-serine
Km_a5_pSer = KM_B0(5);  kcat_a5 = kcat_B0(5);   % 3-P-L-serine-> L-serine (m_in)

%==================
% Block 1: Differentiator
%==================

% Second-substrate constants: the workbook gives one K_M per reaction, so the
% cosubstrate K_M values below are assumed, as for K_M of the cancelers.
Km_a2_asuc = 0.05;  Km_a2_ppi = 0.1;            % argininosuccinate, diphosphate
Km_a3_asp  = 0.3;   Km_a3_gln = 0.3;            % L-aspartate, L-glutamine
Km_a4_php  = 0.01;                              % 3-phospho-hydroxypyruvate
m_asuc_B0 = 0.05; m_ppi_B0 = 0.5;  m_asp_B0 = 1; m_gln_B0 = 2; m_php_B0 = 0.05;
K_EV_AC = 10*vesicles(row,1);            % activation constant in DATASET EV units (assumed)

reactionNames = {'r_Sf','r_Sd','r_Ff','r_Fd','r_Fr','r_I'};

KM = nan(1,5);  kcat = nan(1,5);
colKcat = colOf("k_cat");
for i = 1:5
    r = find(cellfun(@(x) isText(x) && strcmp(strtrim(char(x)),reactionNames{i}),raw(:,2)),1);
    assert(~isempty(r),'Reaction %s not found in column B of the workbook.',reactionNames{i});
    KM(i)   = raw{r,colKM};
    %kcat(i) = raw{r,colSA}*raw{r,colMW}/60;
    kcat(i) = raw{r,colKcat};
end
assert(all(isfinite([KM kcat])) && all([KM kcat] > 0),'K_M, SA and MW of the five reactions must be positive numbers.');

KM_Sf  = KM(1);  kcat_Sf  = kcat(1);       % r_S :  m_in -> m_S
KM_Sd = KM(2);  kcat_Sd = kcat(2);       % r_Sr:  m_S  -> m_d1  (slow rail)
KM_Ff  = KM(3);  kcat_Ff  = kcat(3);       % r_F :  m_in -> m_F
KM_Fd = KM(4);  kcat_Fd = kcat(4);       % r_Fr:  m_F  -> m_P1  (fast rail)
KM_I  = KM(5);  kcat_I  = kcat(5);       % r_I :  m_P1 + m_d1 -> products (K_M used for both substrates: assumed)

%==================
% Block 2: Differentiator
%==================

reactionNames2 = {'r_Sf_2','r_Sd_2','r_Ff_2','r_Fd_2','r_I_2'};
KM2 = nan(1,5);  kcat2 = nan(1,5);
colKcat = colOf("k_cat");
for i = 1:5
    r = find(cellfun(@(x) isText(x) && strcmp(strtrim(char(x)),reactionNames2{i}),raw(:,2)),1);
    assert(~isempty(r),'Reaction %s not found in column B of the workbook.',reactionNames2{i});
    KM2(i)   = raw{r,colKM};
    % kcat2(i) = raw{r,colSA}*raw{r,colMW}/60;
    kcat2(i) = raw{r,colKcat};   
    
end
KM_Sf_2  = KM2(1);  kcat_Sf_2  = kcat2(1);   % r_S_2 : m_d1 -> m_S2
KM_Sd_2 = KM2(2);  kcat_Sd_2 = kcat2(2);   % r_Sr_2: m_S2 -> m_N2 (slow rail)
KM_Ff_2  = KM2(3);  kcat_Ff_2  = kcat2(3);   % r_F_2 : m_d1 -> m_F2
KM_Fd_2 = KM2(4);  kcat_Fd_2 = kcat2(4);   % r_Fr_2: m_F2 -> m_d2 (fast rail)
KM_I_2  = KM2(5);  kcat_I_2  = kcat2(5);   % r_I_2 : m_d2 + m_N2 -> products (K_M for both: assumed)

%==================
% Block 3: Counter
%==================

reactionNamesC = {'r_c','r_I2','r_I1'};
KMC = nan(1,3);  kcatC = nan(1,3);
for i = 1:3
    r = find(cellfun(@(x) isText(x) && strcmp(strtrim(char(x)),reactionNamesC{i}),raw(:,2)),1);
    assert(~isempty(r),'Reaction %s not found in column B of the workbook.',reactionNamesC{i});
    KMC(i)   = raw{r,colKM};
    % kcatC(i) = raw{r,colSA}*raw{r,colMW}/60;
    kcatC(i) = raw{r,colKcat};
end
%assert(all(isfinite([KM2 kcat2 KMC kcatC])) && all([KM2 kcat2 KMC kcatC] > 0),'K_M, SA and MW must be positive numbers.');
KM_Cp   = KMC(1);  kcat_Cp = kcatC(1);   % r_c  = r_Cp: buffered substrate -> m_c, activated by m_d2
KM_Cg_c = KMC(2);  kcat_Cg = kcatC(2);   % r_I1 = r_Cg: m_c + m_I -> products;   K_M for m_c
KM_Cg_I = KMC(2);                        %               K_M for m_I: assumed equal
KM_Cr   = KMC(3);  kcat_Cr = kcatC(3);   % r_I2 = r_Cr: m_I -> products, activated by m_d2

%==================
% Block 3: NN
%==================

% Neuron: the same sheet. Weight r_Ni, removal r_NiDecay; the bias source uses the row r_Nb.
reactionNamesN = {'r_Nw','r_Nr','r_Nb'};
KMN = nan(1,3);  kcatN = nan(1,3);
for i = 1:3
    r = find(cellfun(@(x) isText(x) && strcmp(strtrim(char(x)),reactionNamesN{i}),raw(:,2)),1);
    assert(~isempty(r),'Reaction %s not found in column B of the workbook.',reactionNamesN{i});
    KMN(i)   = raw{r,colKM};
    % kcatN(i) = raw{r,colSA}*raw{r,colMW}/60;
    kcatN(i) = raw{r,colKcat};
end
%assert(all(isfinite([KMN kcatN])) && all([KMN kcatN] > 0),'K_M, SA and MW of the neuron reactions must be positive numbers.');
KM_Nw = KMN(1);  kcat_Nw = kcatN(1);     % r_Nw = r_Ni     : weight (its K_M is used for m_c: assumed)
KM_Nr = KMN(2);  kcat_Nr = kcatN(2);     % r_Nr = r_NiDecay: removal of m_out
KM_Nb = KMN(3);  kcat_Nb = kcatN(3);     % r_Nb            : bias source


%%
%\[text] ## 3. Design choices
E_max = 0.02;                            % largest feasible enzyme concentration [mM] (20 uM)
m_in_target = 0.00125;                   % resting input level [mM]

% Differentiator 1
tau_F_1 = 30;                             % delay of the fast branch (pool m_F) [s]
tau_S_1 = 60;                             % delay of the slow branch (pool m_S) [s];  Delta_t = tau_S1 - tau_F1
tau_r_1 = 30;                             % common removal time of the two rails m_P1 and m_d1 [s]
w_1    = 1;                              % gain: m_P1 - m_d1 = w_1*(fast-delayed m_in - slow-delayed m_in)
% Differentiator 2
tau_F_2 = 30;                             % delay of the fast branch (pool m_F2) [s]
tau_S_2 = 60;                             % delay of the slow branch (pool m_S2) [s]
tau_r_2 = 30;                             % common removal time of the two rails m_d2 and m_N2 [s]
%                                          (its gain is fixed by the load: w_2 = tau_r2/(2*tau_r1))
% Canceler r_I and r_I_2
ideal_canceler = true;                   % true : strong enough to leave the rails at rest_1, rest_2
%                                          false: limited to E_max, to see what is feasible
rest_1 = 1e-7;                           % wanted resting level of the rails of D1 [mM] (ideal canceler)
rest_2 = 1e-8;                           % wanted resting level of the rails of D2 [mM] (ideal canceler)
% Counter
share_c       = 0.9;                     % share of the removal of m_d2 that goes through r_c (-> m_c);
%                                          the rest goes through r_I2. 1 = r_I2 switched off.
tau_reset     = 600;                     % reset time constant of m_c [s]
m_I_base      = 2;                       % buffered level of m_I [mM] (illustrative)
tau_I_supply  = 60;                      % refill time of m_I [s]
% Neuron
w_N   = 1;                               % weight: m_out = w_N*m_c + b_N
b_N   = 5e-7;                            % bias [mM]
tau_N = 30;                              % response time of the neuron [s]
m_b   = 0.01;                            % buffered bias species [mM] (illustrative)
consuming_N = false;                     % false = m_c activates r_Nw; true = r_Nw consumes m_c
assert(w_N >= 0 && b_N >= 0 && tau_N > 0,'The neuron needs w_N >= 0, b_N >= 0 and tau_N > 0.');
assert(tau_S_1 > tau_F_1 && tau_S_2 > tau_F_2,'The slow branch must have the larger delay.');
assert(share_c > 0 && share_c <= 1,'share_c must lie in (0,1].');
%%
%\[text] ## 4. Rate constants, rest state and enzyme concentrations
% ---- Differentiator 1 ----
k_r_1 = 1/tau_r_1;                         % removal of each rail
k_Ff_1  = w_1*k_r_1;                         % m_in -> m_F
k_Sf_1  = k_Ff_1;                              % m_in -> m_S: the same flux into both branches
k_Fd_1 = 1/tau_F_1;                         % m_F -> m_P1
k_Sd_1 = 1/tau_S_1;                         % m_S -> m_d1
k_Fr_1 = k_r_1;                             % m_P1 -> products (new reaction)
k_F_2 = k_r_1/2;  k_S_2 = k_r_1/2;         % m_d1 -> m_F2 and m_d1 -> m_S2: together the k_r of m_d1

flux_1 = k_Ff_1*m_in_target;                % resting flux through each branch [mM/s]
m_F_ref = k_Ff_1*m_in_target/k_Fd_1;
m_S_ref = k_Sf_1*m_in_target/k_Sd_1;
if ideal_canceler
    assert(k_r_1*rest_1 < flux_1,'rest_1 is too high: lower it.');
    k_I = (flux_1-k_r_1*rest_1)/rest_1^2;
else
    k_I = kcat_I*E_max/(KM_I*KM_I);
end
p1_ref = (-k_r_1+sqrt(k_r_1^2+4*k_I*flux_1))/(2*k_I);      % resting m_P1 = resting m_d1

% ---- Differentiator 2 (input: m_d1) ----
k_r_2   = 1/tau_r_2;
k_Fr_2 = 1/tau_F_2;                       % m_F2 -> m_d2
k_Sr_2 = 1/tau_S_2;                       % m_S2 -> m_N2
k_c    = share_c*k_r_2;                   % m_d2 -> m_c            (counter, r_c)
k_I2c  = (1-share_c)*k_r_2;               % m_d2 + m_I -> products (counter, r_I2), at m_I = m_I_base
k_Sx_2 = k_r_2;                           % removal of m_N2 (new reaction)
flux_2 = k_F_2*p1_ref;                   % resting flux through each branch of D2
m_F2_ref = k_F_2*p1_ref/k_Fr_2;
m_S2_ref = k_S_2*p1_ref/k_Sr_2;
if ideal_canceler
    assert(k_r_2*rest_2 < flux_2,'rest_2 is too high for the resting m_d1: lower it.');
    k_I_2 = (flux_2-k_r_2*rest_2)/rest_2^2;
else
    k_I_2 = kcat_I_2*E_max/(KM_I_2*KM_I_2);
end
p2_ref = (-k_r_2+sqrt(k_r_2^2+4*k_I_2*flux_2))/(2*k_I_2);  % resting m_d2 = resting m_N2

% ---- Enzyme concentrations, from k = kcat*[E]/K_M ----
names_all = {'r_F';'r_Fr';'r_S';'r_Sr';'r_I';'r_Fx (new)'; ...
             'r_F_2';'r_Fr_2';'r_S_2';'r_Sr_2';'r_I_2';'r_Sx_2 (new)';'r_c = r_Cp';'r_I2 = r_Cr'};
k_all = [k_Ff_1;k_Fd_1;k_Sf_1;k_Sd_1;k_I;k_Fr_1;k_F_2;k_Fr_2;k_S_2;k_Sr_2;k_I_2;k_Sx_2;k_c;k_I2c];
E_all = [k_Ff_1*KM_Ff/kcat_Ff; k_Fd_1*KM_Fd/kcat_Fd; k_Sf_1*KM_Sf/kcat_Sf; k_Sd_1*KM_Sd/kcat_Sd; k_I*KM_I^2/kcat_I; ...
         k_Fr_1*KM_Fd/kcat_Fd; ...
         k_F_2*KM_Ff_2/kcat_Ff_2; k_Fr_2*KM_Fd_2/kcat_Fd_2; k_S_2*KM_Sf_2/kcat_Sf_2; k_Sr_2*KM_Sd_2/kcat_Sd_2; ...
         k_I_2*KM_I_2^2/kcat_I_2; k_Sx_2*KM_Fd_2/kcat_Fd_2; ...
         k_c*KM_Cp/kcat_Cp; k_I2c*KM_Cr/(kcat_Cr*(m_I_base/(KM_Cr+m_I_base)))];
reaction_parameters = table(names_all,k_all,1e3*E_all, ...
    'VariableNames',{'Reaction','k_per_s_or_per_mM_s','Enzyme_uM'});
disp(reaction_parameters);
for i = find(E_all > E_max).'
    warning('%s needs %.3g uM, above the limit of %.3g uM.',names_all{i},1e3*E_all(i),1e3*E_max);
end
fprintf('Rest: m_F %.3g, m_S %.3g, rails of D1 %.3g, m_F2 %.3g, m_S2 %.3g, rails of D2 %.3g mM\n', ...
    m_F_ref,m_S_ref,p1_ref,m_F2_ref,m_S2_ref,p2_ref);

% ---- Counter (consumes m_d2) ----
k_I_supply = 1/tau_I_supply;
satI_base  = m_I_base/(KM_Cr+m_I_base);                  % r_I2 is sized at m_I = m_I_base
vc_ref  = k_c*p2_ref;                                    % resting production of m_c [mM/s]
I_balance = @(x) k_I_supply*(m_I_base-x)-k_I2c*p2_ref*(x/(KM_Cr+x))/satI_base-vc_ref;
m_I_ref = fzero(I_balance,[0 m_I_base]);
% Neuron rate constants. If r_Nw consumes m_c it is part of the reset of the counter.
k_o = 1/tau_N;   k_w = w_N*k_o;   k_b = b_N*k_o/m_b;
k_reset_I1 = 1/tau_reset-consuming_N*k_w;                % what r_I1 still has to provide
assert(k_reset_I1 > 0,'With consuming_N = true the neuron alone resets m_c faster than tau_reset: lower w_N/tau_N or tau_reset.');
E_Cg = k_reset_I1*KM_Cg_c/(kcat_Cg*(m_I_ref/(KM_Cg_I+m_I_ref)));   V_Cg = kcat_Cg*E_Cg;
c_balance = @(x) vc_ref-V_Cg*x/(KM_Cg_c+x)*m_I_ref/(KM_Cg_I+m_I_ref)-consuming_N*k_w*x;
m_c_ref = fzero(c_balance,[0 vc_ref*tau_reset*10]);
m_out_ref = w_N*m_c_ref+b_N;
E_Nw = k_w*KM_Nw/kcat_Nw;  E_Nb = k_b*KM_Nb/kcat_Nb;  E_Nr = k_o*KM_Nr/kcat_Nr;
fprintf('Neuron: enzymes r_Nw %.3g, r_Nb %.3g, r_Nr %.3g uM;  rest m_out %.3g mM\n',1e3*E_Nw,1e3*E_Nb,1e3*E_Nr,m_out_ref);
fprintf('Counter: r_I1 %.3g uM;  rest m_c %.3g mM, m_I %.6g mM\n',1e3*E_Cg,m_c_ref,m_I_ref);

% Block 0 supplies what r_S and r_F withdraw from m_in.
j_S_ref = k_Sf_1*m_in_target;  j_F_ref = k_Ff_1*m_in_target;
%%
%\[text] ## 5. Block 0: kinetics from BRENDA, design concentrations and capacities
%\[text] The pool levels are no longer free: with K\_M fixed by the enzyme, the resting
%\[text] concentration that gives the wanted response time tau\_B0 follows from the
%\[text] linearised balance d(P)/dt = J - V\*P/(K\_M+P) about the operating point.
tau_B0 = 30;                                     % response time of each Block 0 pool [s]



% ---- Resting pool levels follow from K_M and the wanted response time ----
J_B0_ref = j_S_ref + j_F_ref;                    % L-serine supply = what r_Sf and r_Ff withdraw
poolRef  = @(KM) (-KM + sqrt(KM^2 + 4*J_B0_ref*KM*tau_B0))/2;
cAMP_ref = poolRef(Km_a1_cAMP);                  % each pool is set by its consuming reaction
AMP_ref  = poolRef(Km_a2_AMP);
ATP_ref  = poolRef(Km_a3_ATP);
glu_ref  = poolRef(Km_a4_glu);
pSer_ref = poolRef(Km_a5_pSer);

% Buffered cosubstrates enter as constant saturation fractions.
a2_fraction = (m_asuc_B0/(Km_a2_asuc+m_asuc_B0))*(m_ppi_B0/(Km_a2_ppi+m_ppi_B0));
a3_fraction = (m_asp_B0/(Km_a3_asp+m_asp_B0))*(m_gln_B0/(Km_a3_gln+m_gln_B0));
a4_fraction = m_php_B0/(Km_a4_php+m_php_B0);

% ---- Capacities and enzyme concentrations ----
% r_EV is the receptor-coupled step that turns the EV count into cAMP production.
% It has no row in the workbook, so its capacity is an assumption, not a measurement.
EV_activity_ref = vesicles(row,1)/(K_EV_AC+vesicles(row,1));
V_EV = J_B0_ref/EV_activity_ref;
V_a1 = J_B0_ref/(cAMP_ref/(Km_a1_cAMP+cAMP_ref));
V_a2 = J_B0_ref/((AMP_ref/(Km_a2_AMP+AMP_ref))*a2_fraction);
V_a3 = J_B0_ref/((ATP_ref/(Km_a3_ATP+ATP_ref))*a3_fraction);
V_a4 = J_B0_ref/((glu_ref/(Km_a4_glu+glu_ref))*a4_fraction);
V_a5 = J_B0_ref/(pSer_ref/(Km_a5_pSer+pSer_ref));

V_B0    = [V_a1;V_a2;V_a3;V_a4;V_a5];
E_B0    = V_B0./kcat_B0(:);                                  % [mM]
pool_B0 = [cAMP_ref;AMP_ref;ATP_ref;glu_ref;pSer_ref];       % pool consumed by each reaction

block0_parameters = table(reactionNamesB0(:), KM_B0(:), kcat_B0(:), ...
    1e3*E_B0, V_B0, pool_B0, pool_B0./KM_B0(:), ...
    'VariableNames',{'Reaction','Km_mM','kcat_per_s','Enzyme_uM', ...
                     'Vmax_mM_per_s','Pool_mM','Pool_over_Km'});
disp(block0_parameters);
for i = find(E_B0 > E_max).'
    warning('%s needs %.3g uM, above the limit of %.3g uM.',reactionNamesB0{i},1e3*E_B0(i),1e3*E_max);
end
fprintf('Block 0 total enzyme: %.3g uM;  largest pool / K_M: %.2g;  assumed r_EV capacity %.3g mM/s\n', ...
    1e3*sum(E_B0), max(pool_B0./KM_B0(:)), V_EV);
%%
%[text] 6\. Simulation (analysed distance) from the rest state
% Each step: (1) production and first-order removal with explicit Euler,
% giving trial values P and Q of the two rails; (2) the canceler applied
% implicitly. With D = P - Q (which the canceler does not change), the new
% P solves   P_new = P - dt*k*P_new*(P_new - D),   a quadratic with one
% positive root; then Q_new = P_new - D. Stable for any canceler strength.
EV = vesicles(row,:).';
m_cAMP = cAMP_ref*ones(N,1); m_AMP = AMP_ref*ones(N,1); m_ATP = ATP_ref*ones(N,1);
m_glu = glu_ref*ones(N,1);   m_pSer = pSer_ref*ones(N,1); m_in = m_in_target*ones(N,1);
m_F  = m_F_ref*ones(N,1);    m_S  = m_S_ref*ones(N,1);  m_P1 = p1_ref*ones(N,1);  m_d1 = p1_ref*ones(N,1);
m_F2 = m_F2_ref*ones(N,1);   m_S2 = m_S2_ref*ones(N,1); m_d2 = p2_ref*ones(N,1);  m_N2 = p2_ref*ones(N,1);
m_c = m_c_ref*ones(N,1);     m_I = m_I_ref*ones(N,1);   m_out = m_out_ref*ones(N,1);

for n = 1:N-1
    dt = t(n+1)-t(n);
    % Block 0
    v_EV = V_EV*EV(n)/(K_EV_AC+EV(n));                               % EV count -> cAMP (assumed)
    v_a1 = V_a1*m_cAMP(n)/(Km_a1_cAMP+m_cAMP(n));                    % cAMP         -> AMP
    v_a2 = V_a2*m_AMP(n)/(Km_a2_AMP+m_AMP(n))*a2_fraction;           % AMP          -> ATP
    v_a3 = V_a3*m_ATP(n)/(Km_a3_ATP+m_ATP(n))*a3_fraction;           % ATP          -> L-glutamate
    v_a4 = V_a4*m_glu(n)/(Km_a4_glu+m_glu(n))*a4_fraction;           % L-glutamate  -> 3-P-L-serine
    v_a5 = V_a5*m_pSer(n)/(Km_a5_pSer+m_pSer(n));                    % 3-P-L-serine -> L-serine

    m_cAMP(n+1) = m_cAMP(n)+dt*(v_EV-v_a1);
    m_AMP(n+1)  = m_AMP(n) +dt*(v_a1-v_a2);
    m_ATP(n+1)  = m_ATP(n) +dt*(v_a2-v_a3);
    m_glu(n+1)  = m_glu(n) +dt*(v_a3-v_a4);
    m_pSer(n+1) = m_pSer(n)+dt*(v_a4-v_a5);
    m_in(n+1)   = m_in(n)  +dt*(v_a5-(k_Sf_1+k_Ff_1)*m_in(n));

    % ---- Differentiator 1 ----
    m_F(n+1) = m_F(n)+dt*(k_Ff_1*m_in(n)-k_Fd_1*m_F(n));                    % fast delay pool
    m_S(n+1) = m_S(n)+dt*(k_Sf_1*m_in(n)-k_Sd_1*m_S(n));                    % slow delay pool
    P = m_P1(n)+dt*(k_Fd_1*m_F(n)-k_Fr_1*m_P1(n));                         % fast rail
    Q = m_d1(n)+dt*(k_Sd_1*m_S(n)-(k_F_2+k_S_2)*m_d1(n));                % slow rail; removal = load of D2
    D = P-Q;                                                           % canceler, implicit over dt
    qb = 1-dt*k_I*D;  qd = sqrt(qb^2+4*dt*k_I*P);
    if qb > 0
        P = 2*P/(qb+qd);
    else
        P = (qd-qb)/(2*dt*k_I);
    end
    m_P1(n+1) = max(P,0);  m_d1(n+1) = max(P-D,0);

    % ---- Differentiator 2 (fed by the m_d1 it consumes) ----
    m_F2(n+1) = m_F2(n)+dt*(k_F_2*m_d1(n)-k_Fr_2*m_F2(n));             % fast delay pool
    m_S2(n+1) = m_S2(n)+dt*(k_S_2*m_d1(n)-k_Sr_2*m_S2(n));             % slow delay pool
    kI2_now = k_I2c*(m_I(n)/(KM_Cr+m_I(n)))/satI_base;                 % r_I2 at the present m_I
    P = m_d2(n)+dt*(k_Fr_2*m_F2(n)-(k_c+kI2_now)*m_d2(n));             % fast rail = output; removal = counter
    Q = m_N2(n)+dt*(k_Sr_2*m_S2(n)-k_Sx_2*m_N2(n));                    % slow rail
    D = P-Q;                                                           % canceler, implicit over dt
    qb = 1-dt*k_I_2*D;  qd = sqrt(qb^2+4*dt*k_I_2*P);
    if qb > 0
        P = 2*P/(qb+qd);
    else
        P = (qd-qb)/(2*dt*k_I_2);
    end
    m_d2(n+1) = max(P,0);  m_N2(n+1) = max(P-D,0);

    % ---- Counter (consumes m_d2) ----
    v_c  = k_c*m_d2(n);                                                % r_c : m_d2 -> m_c
    v_I2 = kI2_now*m_d2(n);                                            % r_I2: m_d2 + m_I -> products
    v_Cg = V_Cg*m_c(n)/(KM_Cg_c+m_c(n))*m_I(n)/(KM_Cg_I+m_I(n));       % r_I1: m_c + m_I -> products
    v_Nw = k_w*m_c(n);                                                 % r_Nw: weight
    m_c(n+1) = m_c(n)+dt*(v_c-v_Cg-consuming_N*v_Nw);
    m_I(n+1) = m_I(n)+dt*(k_I_supply*(m_I_base-m_I(n))-v_I2-v_Cg);

    % ---- Neuron ----
    m_out(n+1) = m_out(n)+dt*(v_Nw+k_b*m_b-k_o*m_out(n));              % weight + bias - removal

    pools = [m_out(n+1),m_in(n+1),m_F(n+1),m_S(n+1),m_P1(n+1),m_d1(n+1),m_F2(n+1),m_S2(n+1),m_d2(n+1),m_N2(n+1),m_c(n+1),m_I(n+1)];
    assert(all(isfinite(pools)) && all(pools >= 0),'Invalid step at n = %d; reduce dt0.',n);
end
%%
%\[text] ## 7. Diagnostics
% The identity behind the design, computed by plain filtering of m_in:
% m_P1 - m_d1 = w_1*(m_in delayed by tau_F1 - m_in delayed by tau_S1), smoothed with tau_r1.
x_F = m_in;  x_S = m_in;  ideal_1 = zeros(N,1);
for n = 1:N-1
    x_F(n+1) = x_F(n)+dt0*(m_in(n)-x_F(n))/tau_F_1;
    x_S(n+1) = x_S(n)+dt0*(m_in(n)-x_S(n))/tau_S_1;
    ideal_1(n+1) = ideal_1(n)+dt0*(w_1*(x_F(n)-x_S(n))-ideal_1(n))/tau_r_1;
end
diff_1 = m_P1-m_d1;
fprintf('\nD1: m_P1 - m_d1 against the ideal subtraction: largest gap %.3g mM (signal up to %.3g mM)\n', ...
    max(abs(diff_1-ideal_1)),max(abs(ideal_1)));
fprintf('D1 output m_d1: rest %.3g, peak %.3g at %.3g min, last %.3g mM;  m_P1 peak %.3g mM\n', ...
    m_d1(1),max(m_d1),t(find(m_d1 == max(m_d1),1))/60,m_d1(end),max(m_P1));
fprintf('D2 output m_d2: rest %.3g, peak %.3g at %.3g min, last %.3g mM;  m_N2 peak %.3g mM\n', ...
    m_d2(1),max(m_d2),t(find(m_d2 == max(m_d2),1))/60,m_d2(end),max(m_N2));
fprintf('Counter m_c: rest %.3g, peak %.3g (%.3g times the rest) at %.3g min, last %.3g mM;  m_I lowest %.6g mM\n', ...
    m_c(1),max(m_c),max(m_c)/m_c(1),t(find(m_c == max(m_c),1))/60,m_c(end),min(m_I));
target_N = w_N*m_c+b_N;
fprintf('Neuron m_out: rest %.3g, peak %.3g at %.3g min, last %.3g mM;  largest gap to w_N*m_c + b_N: %.3g mM\n', ...
    m_out(1),max(m_out),t(find(m_out == max(m_out),1))/60,m_out(end),max(abs(m_out-target_N)));
fprintf('Largest pool / K_M (first-order kinetics need << 1): m_F %.2g, m_S %.2g, m_F2 %.2g, m_S2 %.2g\n', ...
    max(m_F)/KM_Fd,max(m_S)/KM_Sd,max(m_F2)/KM_Fd_2,max(m_S2)/KM_Sd_2);
%%
%\[text] ## 8. Figures
% 1. Input.
figure;
plot(t/60,m_in,'LineWidth',2); grid on;
xlabel('Time [min]'); ylabel('L-serine: m_{in} [mM]'); set(gca,'FontSize',fontsize+4);

% 2. Differentiator 1: the two delay pools.
figure;
plot(t/60,[m_F/m_F(1),m_S/m_S(1),m_in/m_in(1)],'LineWidth',2); grid on;
xlabel('Time [min]'); ylabel('relative to rest'); xlim([0,min(30,t(end)/60)]);
legend('m_F (fast delay pool)','m_S (slow delay pool)','m_{in}','Location','best'); set(gca,'FontSize',fontsize);

% 3. Differentiator 1: the two rails, and the falling slope of the input.
figure;
yyaxis left;  plot(t/60,[m_d1,m_P1],'LineWidth',2); ylabel('Concentration [mM]');
yyaxis right; plot(t/60,-gradient(m_in,t),':','LineWidth',1.5); ylabel('-dm_{in}/dt [mM/s]');
xlabel('Time [min]'); xlim([0,min(150,t(end)/60)]); grid on;
legend('m_{d1} (slow rail, output)','m_{P1} (fast rail)','input falling slope','Location','best');
set(gca,'FontSize',fontsize);

% 4. Differentiator 2: the two rails, and the slope of m_d1.
figure;
yyaxis left;  plot(t/60,[m_d2,m_N2],'LineWidth',2); ylabel('Concentration [mM]');
yyaxis right; plot(t/60,gradient(m_d1,t),':','LineWidth',1.5); ylabel('dm_{d1}/dt [mM/s]');
xlabel('Time [min]'); xlim([0,min(150,t(end)/60)]); grid on;
legend('m_{d2} (fast rail, output)','m_{N2} (slow rail)','slope of m_{d1}','Location','best');
set(gca,'FontSize',fontsize);

% 5. Check of the identity: difference of the rails against the ideal subtraction.
figure;
plot(t/60,diff_1,'LineWidth',2); hold on; plot(t/60,ideal_1,'--','LineWidth',1.5); grid on;
xlabel('Time [min]'); ylabel('[mM]'); xlim([0,min(150,t(end)/60)]);
legend('m_{P1} - m_{d1}','w_1 (fast-delayed - slow-delayed m_{in})','Location','best'); set(gca,'FontSize',fontsize);

% 6. Counter.
figure;
yyaxis left;  plot(t/60,m_c,'LineWidth',2); ylabel('m_c [mM]');
yyaxis right; plot(t/60,m_d2,'--','LineWidth',1.5); ylabel('m_{d2} [mM]');
xlabel('Time [min]'); xlim([0,min(150,t(end)/60)]); grid on;
legend('m_c (counter)','m_{d2} (consumed by the counter)','Location','best'); set(gca,'FontSize',fontsize);

% 7. Neuron: output against w_N*m_c + b_N.
figure;
plot(t/60,m_out,'LineWidth',2); hold on; plot(t/60,target_N,'--','LineWidth',1.5);
plot(t/60,b_N*ones(N,1),':','LineWidth',1.2); grid on;
xlabel('Time [min]'); ylabel('[mM]'); xlim([0,min(150,t(end)/60)]);
legend('m_{out} (neuron)','w_N m_c + b_N','bias b_N','Location','best'); set(gca,'FontSize',fontsize);
%%
%\[text] ## Path helper (unchanged from your original script)
function addingPathParentFolderByName(targetName)
    % Add the parent folder called targetName (and its subfolders) to the path.
    currFolder = pwd;
    found = false;
    while true
        [parentFolder, currentName] = fileparts(currFolder);
        if strcmpi(currentName, targetName)
            found = true;
            break;
        end
        if isempty(parentFolder) || strcmp(currFolder, parentFolder)
            break;
        end
        currFolder = parentFolder;
    end
    if found
        addpath(genpath(currFolder));
        fprintf('Adding matlab path to: %s\n', currFolder);
    else
        error('Folder named "%s" not found in any parent directory.', targetName);
    end
end

%\[text] ## Kinetics helper
function [KM, kcat] = readKinetics(raw, names, colKM, colKcat)
    % Read K_M [mM] and k_cat [1/s] for the listed reaction labels from column B.
    isText = @(x) ischar(x) || (isstring(x) && isscalar(x) && ~ismissing(x));
    n = numel(names);  KM = nan(1,n);  kcat = nan(1,n);
    for i = 1:n
        r = find(cellfun(@(x) isText(x) && strcmp(strtrim(char(x)),names{i}), raw(:,2)),1);
        assert(~isempty(r),'Reaction %s not found in column B of the workbook.',names{i});
        KM(i)   = raw{r,colKM};
        kcat(i) = raw{r,colKcat};
    end
    assert(all(isfinite([KM kcat])) && all([KM kcat] > 0), ...
        'K_M and k_cat must be positive and finite for: %s', strjoin(names,', '));
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":36.3}
%---
%[output:4c4d511c]
%   data: {"dataType":"text","outputData":{"text":"Adding matlab path to: D:\\code\\bioelectrodynamics\\code\n","truncated":false}}
%---
