clc;
clearvars;
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
%[text] Pulling Enzyme constants from the BRENDA workbook
%[text] Description: This block uses your credentials at "Brenda\\\_credentials.txt" for searching at \[Brenda website\](https://www.brenda-enzymes.info/functional\_parameter\_statistics.php) the enzyme parameters. Your credentials should be included in a txt file with two lines, the first line includes your email, and the second line your password. The enzymes are identified from the file "Block\\\_7\\\_Selected\\\_Motifs.txt", which is completed by the second step and includes the chemical reaction network of the algorithm. After the code reads the enzyme parameters it prepares an excel table with the listed values per enzyme.
% % Uncomment below for reading the constants from https://www.brenda-enzymes.info/functional_parameter_statistics.php 
% 
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
% 
% brendaKineticsTable('Selected_Full_Chain_1.txt','kinetics_BRENDA.xlsx', ...
%         'BrendaFlatFile', which('brenda_2026_1.txt'), 'KcatSource','BRENDA else SAxMW');
%%

% % Sheet "values" of kinetics_BRENDA.xlsx: reaction label in column B, and the
% % columns K_M [mM], SA [umol/min/mg] and MW [mg/umol] found by their titles.
% % kcat [1/s] = SA*MW/60.
% raw = readcell('kinetics_BRENDA.xlsx','Sheet','values');
% isText = @(x) ischar(x) || (isstring(x) && isscalar(x) && ~ismissing(x));
% hdrRow = find(any(cellfun(@(x) isText(x) && startsWith(string(x),"K_M"),raw),2),1);
% assert(~isempty(hdrRow),'Header row (a cell starting with "K_M") not found in the workbook.');
% colOf = @(title) find(cellfun(@(x) isText(x) && startsWith(string(x),title),raw(hdrRow,:)),1);
% colKM = colOf("K_M");  colSA = colOf("SA");  colMW = colOf("MW"); colKcat = colOf("k_cat");
% 
% %==================
% % Block 0
% %==================
% 
% % ---- Constants from the workbook ----
% % Labels and chain exactly as the sheets "values" and "Reaction mapping" list them:
% %   r_adapter_1  3',5'-cyclic AMP            -> AMP                  (K_M for cAMP)
% %   r_adapter_2  AMP + argininosuccinate+PPi -> ATP                  (K_M for ATP)
% %   r_adapter_3  ATP + aspartate + glutamine -> L-glutamate (m_d1)   (K_M for ATP)
% %   r_adapter_4  L-glutamate + 3-P-hydroxypyruvate -> 3-P-L-serine   (K_M for glutamate)
% %   r_adapter_5  3-P-L-serine                -> L-serine (m_in)      (K_M for 3-P-serine)
% reactionNamesB0 = {'r_CAMP_1'};
% KM_B0 = nan(1,length(reactionNamesB0));
% kcat_B0 = nan(1,length(reactionNamesB0));
% 
% for i = 1:length(reactionNamesB0)
%     r = find(cellfun(@(x) isText(x) && strcmp(strtrim(char(x)),reactionNamesB0{i}),raw(:,2)),1);
%     assert(~isempty(r),'Reaction %s not found in column B of the workbook.',reactionNamesB0{i});
%     KM_B0(i)   = raw{r,colKM};
%     %kcat(i) = raw{r,colSA}*raw{r,colMW}/60;
%     val = raw{r,colKcat};
%     if isempty(val), val = NaN; end
%     kcat_B0(i) = val;
% 
% end
% % assert(all(isfinite([KM_B0 kcat_B0])) && all([KM_B0 kcat_B0] > 0),'K_M, SA and MW of the five reactions must be positive numbers.');
% 
% 
% Km_a1_cAMP = KM_B0(1);  kcat_a1 = kcat_B0(1);   % cAMP        -> AMP
% % Km_a2_AMP  = KM_B0(2);  kcat_a2 = kcat_B0(2);   % AMP         -> ATP
% % Km_a3_ATP  = KM_B0(3);  kcat_a3 = kcat_B0(3);   % ATP         -> L-glutamate
% % Km_a4_glu  = KM_B0(4);  kcat_a4 = kcat_B0(4);   % L-glutamate -> 3-P-L-serine
% % Km_a5_pSer = KM_B0(5);  kcat_a5 = kcat_B0(5);   % 3-P-L-serine-> L-serine (m_in)
% 
% %==================
% % Block 1: Differentiator
% %==================
% 
% % Second-substrate constants: the workbook gives one K_M per reaction, so the
% % cosubstrate K_M values below are assumed, as for K_M of the cancelers.
% Km_a2_asuc = 0.05;  Km_a2_ppi = 0.1;            % argininosuccinate, diphosphate
% Km_a3_asp  = 0.3;   Km_a3_gln = 0.3;            % L-aspartate, L-glutamine
% Km_a4_php  = 0.01;                              % 3-phospho-hydroxypyruvate
% m_asuc_B0 = 0.05; m_ppi_B0 = 0.5;  m_asp_B0 = 1; m_gln_B0 = 2; m_php_B0 = 0.05;
% K_EV_AC = 10*vesicles(row,1);            % activation constant in DATASET EV units (assumed)
% 
% reactionNames = {'r_Sf','r_Sd','r_Ff','r_Fd','r_Fr','r_I'};
% 
% KM = nan(1,length(reactionNames));
% kcat = nan(1,length(reactionNames));
% colKcat = colOf("k_cat");
% for i = 1:length(reactionNames)
%     r = find(cellfun(@(x) isText(x) && strcmp(strtrim(char(x)),reactionNames{i}),raw(:,2)),1);
%     assert(~isempty(r),'Reaction %s not found in column B of the workbook.',reactionNames{i});
%     KM(i)   = raw{r,colKM};
%     %kcat(i) = raw{r,colSA}*raw{r,colMW}/60;
%     kcat(i) = raw{r,colKcat};
% end
% assert(all(isfinite([KM kcat])) && all([KM kcat] > 0),'K_M, SA and MW of the five reactions must be positive numbers.');
% 
% KM_Sf_1  = KM(1);  kcat_Sf_1  = kcat(1);       % r_S :  m_in -> m_S
% KM_Sd_1 = KM(2);  kcat_Sd_1 = kcat(2);       % r_Sr:  m_S  -> m_d1  (slow rail)
% KM_Ff_1  = KM(3);  kcat_Ff_1  = kcat(3);       % r_F :  m_in -> m_F
% KM_Fd_1 = KM(4);  kcat_Fd_1 = kcat(4);       % r_Fr:  m_F  -> m_P1  (fast rail)
% KM_Fr_1 = KM(5); kcat_Fr_1 = kcat(5);
% KM_I_1  = KM(6);  kcat_I_1  = kcat(6);       % r_I :  m_P1 + m_d1 -> products (K_M used for both substrates: assumed)
% 
% 
% %==================
% % Block 2: Differentiator
% %==================
% 
% reactionNames2 = {'r_Sf_2','r_Sd_2','r_Sr_2','r_Ff_2','r_Fd_2','r_I_2'};
% KM2 = nan(1,length(reactionNames2)); 
% kcat2 = nan(1,length(reactionNames2));
% colKcat = colOf("k_cat");
% for i = 1:length(reactionNames2)
%     r = find(cellfun(@(x) isText(x) && strcmp(strtrim(char(x)),reactionNames2{i}),raw(:,2)),1);
%     assert(~isempty(r),'Reaction %s not found in column B of the workbook.',reactionNames2{i});
%     KM2(i)   = raw{r,colKM};
%     % kcat2(i) = raw{r,colSA}*raw{r,colMW}/60;
%     kcat2(i) = raw{r,colKcat};   
% 
% end
% KM_Sf_2  = KM2(1);  kcat_Sf_2  = kcat2(1);   % r_S_2 : m_d1 -> m_S2
% KM_Sd_2 = KM2(2);  kcat_Sd_2 = kcat2(2);   % r_Sr_2: m_S2 -> m_N2 (slow rail)
% KM_Sr_2 = KM2(3);  kcat_Sr_2 = kcat2(3);   % r_Sr_2: m_S2 -> m_N2 (slow rail)
% KM_Ff_2  = KM2(4);  kcat_Ff_2  = kcat2(4);   % r_F_2 : m_d1 -> m_F2
% KM_Fd_2 = KM2(5);  kcat_Fd_2 = kcat2(5);   % r_Fr_2: m_F2 -> m_d2 (fast rail)
% KM_I_2  = KM2(6);  kcat_I_2  = kcat2(6);   % r_I_2 : m_d2 + m_N2 -> products (K_M for both: assumed)
% 
% 
% %==================
% % Block 3: Counter
% %==================
% 
% reactionNamesC = {'r_c','r_I2','r_I1'};
% KMC = nan(1,3);  kcatC = nan(1,3);
% for i = 1:3
%     r = find(cellfun(@(x) isText(x) && strcmp(strtrim(char(x)),reactionNamesC{i}),raw(:,2)),1);
%     assert(~isempty(r),'Reaction %s not found in column B of the workbook.',reactionNamesC{i});
%     KMC(i)   = raw{r,colKM};
%     % kcatC(i) = raw{r,colSA}*raw{r,colMW}/60;
%     kcatC(i) = raw{r,colKcat};
% end
% %assert(all(isfinite([KM2 kcat2 KMC kcatC])) && all([KM2 kcat2 KMC kcatC] > 0),'K_M, SA and MW must be positive numbers.');
% KM_Cp   = KMC(1);  kcat_Cp = kcatC(1);   % r_c  = r_Cp: buffered substrate -> m_c, activated by m_d2
% KM_Cg_c = KMC(2);  kcat_Cg = kcatC(2);   % r_I1 = r_Cg: m_c + m_I -> products;   K_M for m_c
% KM_Cg_I = KMC(2);                        %               K_M for m_I: assumed equal
% KM_Cr   = KMC(3);  kcat_Cr = kcatC(3);   % r_I2 = r_Cr: m_I -> products, activated by m_d2
% 
% %==================
% % Block 3: NN
% %==================
% 
% % Neuron: the same sheet. Weight r_Ni, removal r_NiDecay; the bias source uses the row r_Nb.
% reactionNamesN = {'r_Nw','r_Nr','r_Nb'};
% KMN = nan(1,3);  kcatN = nan(1,3);
% for i = 1:3
%     r = find(cellfun(@(x) isText(x) && strcmp(strtrim(char(x)),reactionNamesN{i}),raw(:,2)),1);
%     assert(~isempty(r),'Reaction %s not found in column B of the workbook.',reactionNamesN{i});
%     KMN(i)   = raw{r,colKM};
%     % kcatN(i) = raw{r,colSA}*raw{r,colMW}/60;
%     kcatN(i) = raw{r,colKcat};
% end
% %assert(all(isfinite([KMN kcatN])) && all([KMN kcatN] > 0),'K_M, SA and MW of the neuron reactions must be positive numbers.');
% KM_Nw = KMN(1);  kcat_Nw = kcatN(1);     % r_Nw = r_Ni     : weight (its K_M is used for m_c: assumed)
% KM_Nr = KMN(2);  kcat_Nr = kcatN(2);     % r_Nr = r_NiDecay: removal of m_out
% KM_Nb = KMN(3);  kcat_Nb = kcatN(3);     % r_Nb            : bias source
%[text] 3\. Design choices
% E_max = 0.02;                            % largest feasible enzyme concentration [mM] (20 uM)
% m_in_target = 0.00125;                   % resting input level [mM]
% 
% % Differentiator 1
% tau_F_1 = 30;                             % delay of the fast branch (pool m_F) [s]
% tau_S_1 = 60;                             % delay of the slow branch (pool m_S) [s];  Delta_t = tau_S1 - tau_F1
% tau_r_1 = 30;                             % common removal time of the two rails m_P1 and m_d1 [s]
% w_1    = 1;                              % gain: m_P1 - m_d1 = w_1*(fast-delayed m_in - slow-delayed m_in)
% % Differentiator 2
% tau_F_2 = 30;                             % delay of the fast branch (pool m_F2) [s]
% tau_S_2 = 60;                             % delay of the slow branch (pool m_S2) [s]
% tau_r_2 = 30;                             % common removal time of the two rails m_d2 and m_N2 [s]
% %                                          (its gain is fixed by the load: w_2 = tau_r2/(2*tau_r1))
% % Canceler r_I and r_I_2
% ideal_canceler = false;                   % true : strong enough to leave the rails at rest_1, rest_2
% %                                          false: limited to E_max, to see what is feasible
% rest_1 = 1e-7;                           % wanted resting level of the rails of D1 [mM] (ideal canceler)
% rest_2 = 1e-8;                           % wanted resting level of the rails of D2 [mM] (ideal canceler)
% % Counter
% share_c       = 0.9;                     % share of the removal of m_d2 that goes through r_c (-> m_c);
% %                                          the rest goes through r_I2. 1 = r_I2 switched off.
% tau_reset     = 600;                     % reset time constant of m_c [s]
% m_I_base      = 2;                       % buffered level of m_I [mM] (illustrative)
% tau_I_supply  = 60;                      % refill time of m_I [s]
% % Neuron
% w_N   = 1;                               % weight: m_out = w_N*m_c + b_N
% b_N   = 5e-7;                            % bias [mM]
% tau_N = 30;                              % response time of the neuron [s]
% m_b   = 0.01;                            % buffered bias species [mM] (illustrative)
% consuming_N = false;                     % false = m_c activates r_Nw; true = r_Nw consumes m_c
% assert(w_N >= 0 && b_N >= 0 && tau_N > 0,'The neuron needs w_N >= 0, b_N >= 0 and tau_N > 0.');
% assert(tau_S_1 > tau_F_1 && tau_S_2 > tau_F_2,'The slow branch must have the larger delay.');
% assert(share_c > 0 && share_c <= 1,'share_c must lie in (0,1].');
%[text] #### Enzyme concentrations
% % ---- Differentiator 1 ----
% k_r_1 = 1/tau_r_1;                         % removal of each rail
% k_Ff_1  = w_1*k_r_1;                         % m_in -> m_F
% k_Sf_1  = k_Ff_1;                              % m_in -> m_S: the same flux into both branches
% k_Fd_1 = 1/tau_F_1;                         % m_F -> m_P1
% k_Sd_1 = 1/tau_S_1;                         % m_S -> m_d1
% k_Fr_1 = k_r_1;                             % m_P1 -> products (new reaction)
% k_Ff_2 = k_r_1/2;
% k_Sf_2 = k_r_1/2;         % m_d1 -> m_F2 and m_d1 -> m_S2: together the k_r of m_d1
% 
% flux_1 = k_Ff_1*m_in_target;                % resting flux through each branch [mM/s]
% m_Ff_ref = k_Ff_1*m_in_target/k_Fd_1;
% m_Sf_ref = k_Sf_1*m_in_target/k_Sd_1;
% 
% if ideal_canceler
%     assert(k_r_1*rest_1 < flux_1,'rest_1 is too high: lower it.');
%     k_I_1 = (flux_1-k_r_1*rest_1)/rest_1^2;
% else
%     k_I_1 = kcat_I_1*E_max/(KM_I_1*KM_I_1);
% end
% m_Fd_ref = (-k_r_1+sqrt(k_r_1^2+4*k_I_1*flux_1))/(2*k_I_1);      % resting m_P1 = resting m_d1
% m_Sd_ref = (-k_r_1+sqrt(k_r_1^2+4*k_I_1*flux_1))/(2*k_I_1);      % resting m_P1 = resting m_d1
% 
% % ---- Differentiator 2 (input: m_d1) ----
% k_r_2   = 1/tau_r_2;
% k_Fd_2 = 1/tau_F_2;                       % m_F2 -> m_d2
% k_Sd_2 = 1/tau_S_2;                       % m_S2 -> m_N2
% k_c    = share_c*k_r_2;                   % m_d2 -> m_c            (counter, r_c)
% k_I2c  = (1-share_c)*k_r_2;               % m_d2 + m_I -> products (counter, r_I2), at m_I = m_I_base
% k_Sr_2 = k_r_2;                           % removal of m_N2 (new reaction)
% flux_2 = k_Ff_2*m_Fd_ref;                   % resting flux through each branch of D2
% m_F2_ref = k_Ff_2*m_Fd_ref/k_Fd_2;
% m_S2_ref = k_Sf_2*m_Fd_ref/k_Sd_2;
% if ideal_canceler
%     assert(k_r_2*rest_2 < flux_2,'rest_2 is too high for the resting m_d1: lower it.');
%     k_I_2 = (flux_2-k_r_2*rest_2)/rest_2^2;
% else
%     k_I_2 = kcat_I_2*E_max/(KM_I_2*KM_I_2);
% end
% p2_ref = (-k_r_2+sqrt(k_r_2^2+4*k_I_2*flux_2))/(2*k_I_2);  % resting m_d2 = resting m_N2
% 
% % ---- Enzyme concentrations, from k = kcat*[E]/K_M ----
% names_all = {'r_Ff_1';'r_Fd_1';'r_Fr_1';'r_Sf_1';'r_Sd_1';'r_I_1'; ... %differentiator 1
%              'r_Ff_2';'r_Fd_2';'r_Sf_2';'r_Sd_2';'r_Sr_2';'r_I_2';... %differentiator 2
%              'r_c = r_Cp';'r_I2 = r_Cr'}; %counter
% k_all = [k_Ff_1;k_Fd_1;k_Fr_1;k_Sf_1;k_Sd_1;k_I_1;...%differentiator 1
%          k_Ff_2;k_Fd_2;k_I_2;k_Sf_2;k_Sd_2;k_Sr_2;...%differentiator 2
%          k_c;k_I2c]; %counter
% 
% E_all = [... %differentiator 1
%          k_Ff_1*KM_Ff_1/kcat_Ff_1; k_Fd_1*KM_Fd_1/kcat_Fd_1; k_Fr_1*KM_Fr_1/kcat_Fr_1;...%fast branch
%          k_Sf_1*KM_Sf_1/kcat_Sf_1; k_Sd_1*KM_Sd_1/kcat_Sd_1; k_I_1*KM_I_1^2/kcat_I_1; ...%slow branch
%          ... %differentiator 2 
%          k_Ff_2*KM_Ff_2/kcat_Ff_2; k_Fd_2*KM_Fd_2/kcat_Fd_2; k_I_2*KM_I_2^2/kcat_I_2; ...%fast branch
%          k_Sf_2*KM_Sf_2/kcat_Sf_2; k_Sd_2*KM_Sd_2/kcat_Sd_2; k_Sr_2*KM_Sr_2/kcat_Sr_2; ...%slow branch 
%          k_c*KM_Cp/kcat_Cp; k_I2c*KM_Cr/(kcat_Cr*(m_I_base/(KM_Cr+m_I_base)))]; %counter
% 
% reaction_parameters = table(names_all,k_all,1e3*E_all, ...
%     'VariableNames',{'Reaction','k_per_s_or_per_mM_s','Enzyme_uM'});
% disp(reaction_parameters);
% for i = find(E_all > E_max).'
%     warning('%s needs %.3g uM, above the limit of %.3g uM.',names_all{i},1e3*E_all(i),1e3*E_max);
% end
% fprintf('Rest: m_F %.3g, m_S %.3g, rails of D1 %.3g, m_F2 %.3g, m_S2 %.3g, rails of D2 %.3g mM\n', ...
%     m_Ff_ref,m_Sf_ref,m_Fd_ref,m_F2_ref,m_S2_ref,p2_ref);
% 
% % ---- Counter (consumes m_d2) ----
% k_I_supply = 1/tau_I_supply;
% satI_base  = m_I_base/(KM_Cr+m_I_base);                  % r_I2 is sized at m_I = m_I_base
% vc_ref  = k_c*p2_ref;                                    % resting production of m_c [mM/s]
% I_balance = @(x) k_I_supply*(m_I_base-x)-k_I2c*p2_ref*(x/(KM_Cr+x))/satI_base-vc_ref;
% m_I_ref = fzero(I_balance,[0 m_I_base]);
% % Neuron rate constants. If r_Nw consumes m_c it is part of the reset of the counter.
% k_o = 1/tau_N;   k_w = w_N*k_o;   k_b = b_N*k_o/m_b;
% k_reset_I1 = 1/tau_reset-consuming_N*k_w;                % what r_I1 still has to provide
% assert(k_reset_I1 > 0,'With consuming_N = true the neuron alone resets m_c faster than tau_reset: lower w_N/tau_N or tau_reset.');
% E_Cg = k_reset_I1*KM_Cg_c/(kcat_Cg*(m_I_ref/(KM_Cg_I+m_I_ref)));   V_Cg = kcat_Cg*E_Cg;
% c_balance = @(x) vc_ref-V_Cg*x/(KM_Cg_c+x)*m_I_ref/(KM_Cg_I+m_I_ref)-consuming_N*k_w*x;
% m_c_ref = fzero(c_balance,[0 vc_ref*tau_reset*10]);
% m_out_ref = w_N*m_c_ref+b_N;
% E_Nw = k_w*KM_Nw/kcat_Nw;  E_Nb = k_b*KM_Nb/kcat_Nb;  E_Nr = k_o*KM_Nr/kcat_Nr;
% fprintf('Neuron: enzymes r_Nw %.3g, r_Nb %.3g, r_Nr %.3g uM;  rest m_out %.3g mM\n',1e3*E_Nw,1e3*E_Nb,1e3*E_Nr,m_out_ref);
% fprintf('Counter: r_I1 %.3g uM;  rest m_c %.3g mM, m_I %.6g mM\n',1e3*E_Cg,m_c_ref,m_I_ref);
% 
% % Block 0 supplies what r_S and r_F withdraw from m_in.
% j_S_ref = k_Sf_1*m_in_target;  j_F_ref = k_Ff_1*m_in_target;
%[text] Block 0: kinetics from BRENDA, design concentrations and capacities
% %[text] The pool levels are no longer free: with K\_M fixed by the enzyme, the resting
% %[text] concentration that gives the wanted response time tau\_B0 follows from the
% %[text] linearised balance d(P)/dt = J - V\*P/(K\_M+P) about the operating point.
% tau_B0 = 30;                                     % response time of each Block 0 pool [s]
% 
% 
% 
% % ---- Resting pool levels follow from K_M and the wanted response time ----
% J_B0_ref = j_S_ref + j_F_ref;                    % L-serine supply = what r_Sf and r_Ff withdraw
% poolRef  = @(KM) (-KM + sqrt(KM^2 + 4*J_B0_ref*KM*tau_B0))/2;
% cAMP_ref = poolRef(Km_a1_cAMP);                  % each pool is set by its consuming reaction
% % AMP_ref  = poolRef(Km_a2_AMP);
% % ATP_ref  = poolRef(Km_a3_ATP);
% % glu_ref  = poolRef(Km_a4_glu);
% % pSer_ref = poolRef(Km_a5_pSer);
% 
% % Buffered cosubstrates enter as constant saturation fractions.
% a2_fraction = (m_asuc_B0/(Km_a2_asuc+m_asuc_B0))*(m_ppi_B0/(Km_a2_ppi+m_ppi_B0));
% a3_fraction = (m_asp_B0/(Km_a3_asp+m_asp_B0))*(m_gln_B0/(Km_a3_gln+m_gln_B0));
% a4_fraction = m_php_B0/(Km_a4_php+m_php_B0);
% 
% % ---- Capacities and enzyme concentrations ----
% % r_EV is the receptor-coupled step that turns the EV count into cAMP production.
% % It has no row in the workbook, so its capacity is an assumption, not a measurement.
% EV_activity_ref = vesicles(row,1)/(K_EV_AC+vesicles(row,1));
% V_EV = J_B0_ref/EV_activity_ref;
% V_a1 = J_B0_ref/(cAMP_ref/(Km_a1_cAMP+cAMP_ref));
% % V_a2 = J_B0_ref/((AMP_ref/(Km_a2_AMP+AMP_ref))*a2_fraction);
% % V_a3 = J_B0_ref/((ATP_ref/(Km_a3_ATP+ATP_ref))*a3_fraction);
% % V_a4 = J_B0_ref/((glu_ref/(Km_a4_glu+glu_ref))*a4_fraction);
% % V_a5 = J_B0_ref/(pSer_ref/(Km_a5_pSer+pSer_ref));
% 
% % V_B0    = [V_a1;V_a2;V_a3;V_a4;V_a5];
% V_B0    = [V_a1];
% E_B0    = V_B0./kcat_B0(:);                                  % [mM]
% % pool_B0 = [cAMP_ref;AMP_ref;ATP_ref;glu_ref;pSer_ref];       % pool consumed by each reaction
% pool_B0 = [cAMP_ref];       % pool consumed by each reaction
% 
% block0_parameters = table(reactionNamesB0(:), KM_B0(:), kcat_B0(:), ...
%     1e3*E_B0, V_B0, pool_B0, pool_B0./KM_B0(:), ...
%     'VariableNames',{'Reaction','Km_mM','kcat_per_s','Enzyme_uM', ...
%                      'Vmax_mM_per_s','Pool_mM','Pool_over_Km'});
% disp(block0_parameters);
% for i = find(E_B0 > E_max).'
%     warning('%s needs %.3g uM, above the limit of %.3g uM.',reactionNamesB0{i},1e3*E_B0(i),1e3*E_max);
% end
% fprintf('Block 0 total enzyme: %.3g uM;  largest pool / K_M: %.2g;  assumed r_EV capacity %.3g mM/s\n', ...
%     1e3*sum(E_B0), max(pool_B0./KM_B0(:)), V_EV);

%[text] Simulation (analysed distance) from the rest state
% Each step: (1) production and first-order removal with explicit Euler,
% giving trial values P and Q of the two rails; (2) the canceler applied
% implicitly. With D = P - Q (which the canceler does not change), the new
% P solves   P_new = P - dt*k*P_new*(P_new - D),   a quadratic with one
% positive root; then Q_new = P_new - D. Stable for any canceler strength.
EV = vesicles(row,:).';
k_EV = 0.02; %assumed from previous catalitic reactions
K_EV_AC = 0.1545;
% m_cAMP = cAMP_ref*ones(N,1);
%m_AMP = AMP_ref*ones(N,1);
% m_ATP = ATP_ref*ones(N,1);
% m_glu = glu_ref*ones(N,1);   m_pSer = pSer_ref*ones(N,1);
% m_in = m_in_target*ones(N,1);
% m_Ff_1  = m_Ff_ref*ones(N,1);    m_Sf_1  = m_Sf_ref*ones(N,1); m_Sd_1  = m_Sd_ref*ones(N,1);  m_Fd_1 = m_Fd_ref*ones(N,1);  m_d1 = m_Fd_ref*ones(N,1);
% m_F2 = m_F2_ref*ones(N,1);   m_S2 = m_S2_ref*ones(N,1); m_d2 = p2_ref*ones(N,1);  m_N2 = p2_ref*ones(N,1);
% m_c = m_c_ref*ones(N,1);     m_I = m_I_ref*ones(N,1);   m_out = m_out_ref*ones(N,1);

%[text] Evaluating gains and delays per reaction
% Target fast-branch concentration gain
w_F = 1;

% Fast-branch relaxation times [s]
tau_Fd = 60;
tau_Fr = 30;

% Desired extra lag of slow production relative to fast enzyme [s]
Delta_t = 30;

% Fast branch
k_Fd = 1 / tau_Fd;
k_Fr = 1 / tau_Fr;
k_Ff = w_F * k_Fr;

% Effective lag of the fast enzyme concentration
T_F = tau_Fd + tau_Fr;

% Slow production must have lag T_F + Delta_t
k_Sd = 1 / (T_F + Delta_t);

% Independent output removal
tau_Ik = 30;                  % Illustrative choice [s]
k_Ik   = 1 / tau_Ik;

% Nominal slow concentration gain WITHOUT catalytic removal
w_S_uncoupled = 1;
k_Sf = w_S_uncoupled * k_Ik;

%Canceller
m_cAMP_rest=0.5;
m_d1_rest =1; 
w_F = k_Ff / k_Fr;
m_Fd_rest = w_F * m_cAMP_rest;
k_I = (k_Sf*m_cAMP_rest/m_d1_rest - k_Ik) / m_Fd_rest;

J_EV_rest = (15*k_EV) * EV(1) / (K_EV_AC + EV(1));
m_cAMP_ref = J_EV_rest / (k_Ff + k_Sf);

m_cAMP = zeros(N,1);
m_cAMP(1) = m_cAMP_ref;
m_Ff = zeros(N,1);
m_Fd = zeros(N,1);
m_Sf = zeros(N,1);
m_d1 = zeros(N,1);

for n = 1:N-1
    % Time interval
    dt_n = t(n+1) - t(n);

    % ---------------------------------------------------------
    % ---- BLOCK 0 ----
    % ---------------------------------------------------------
    v_EV = k_EV*EV(n)/(K_EV_AC+EV(n));                               % EV count -> cAMP (assumed)
    %v_a1 = V_a1*m_cAMP(n)/(Km_a1_cAMP+m_cAMP(n));                    % cAMP         -> AMP


    % ---------------------------------------------------------
    % ---- BLOCK 1 ---- 
    % Reaction fluxes: evaluate all using concentrations at n
    % ---------------------------------------------------------

    % Fast branch
    v_Ff = k_Ff * m_cAMP(n);    % Input -> m_Ff
    v_Fd = k_Fd * m_Ff(n);      % m_Ff -> m_Fd
    v_Fr = k_Fr * m_Fd(n);      % Independent turnover of m_Fd

    % Slow branch
    v_Sf = k_Sf * m_cAMP(n);    % Input -> m_Sf
    v_Sd = k_Sd * m_Sf(n);      % m_Sf -> m_d1
    v_Ik = k_Ik * m_d1(n);      % Independent removal of m_d1

    % Catalytic coupling: m_Fd accelerates removal of m_d1
    v_I = k_I * m_Fd(n) * m_d1(n);

    % ---------------------------------------------------------
    % Concentration updates: explicit Euler
    % ---------------------------------------------------------

    % ---------------------------------------------------------
    % ---- BLOCK 0 ----
    % ---------------------------------------------------------

    m_cAMP(n+1) = m_cAMP(n) + dt_n*(v_EV - v_Ff);

    % ---------------------------------------------------------
    % ---- BLOCK 1 ----
    % ---------------------------------------------------------

    m_Ff(n+1) = m_Ff(n) + dt_n * (v_Ff - v_Fd);

    % r_I does not consume its catalyst m_Fd
    m_Fd(n+1) = m_Fd(n) + dt_n * (v_Fd - v_Fr);

    m_Sf(n+1) = m_Sf(n) + dt_n * (v_Sf - v_Sd);

    % Output of the differentiator
    m_d1(n+1) = m_d1(n) + dt_n * (v_Sd - v_Ik - v_I);
end

%[text] Figures
% 1. Input.
figure; %[output:96547eb7]
plot(t/60,EV,'LineWidth',2); grid on; %[output:96547eb7]
ylabel('EV'); %[output:96547eb7]
yyaxis right %[output:96547eb7]
plot(t/60,m_cAMP,'LineWidth',2); %[output:96547eb7]
ylabel('cAMP'); %[output:96547eb7]
xlabel('Time [min]');  %[output:96547eb7]
set(gca,'FontSize',fontsize+4); %[output:96547eb7]

% 2. Differentiator 1: the two delay pools.
figure; %[output:78729cd3]
plot(t/60,m_cAMP,'LineWidth',2); grid on; hold on; %[output:78729cd3]
plot(t/60,m_d1,'LineWidth',2); %[output:78729cd3]
xlabel('Time [min]','Interpreter','latex'); %[output:78729cd3]
ylabel('[mM]','Interpreter','latex'); %[output:78729cd3]
xlim([0,min(30,t(end)/60)]); %[output:78729cd3]
legend('$m_{in}$','$m_{d1}$','Location','best','Interpreter','latex'); %[output:78729cd3]
set(gca,'FontSize',fontsize); %[output:78729cd3]

% % 3. Differentiator 1: the two rails, and the falling slope of the input.
% figure;
% yyaxis left;  plot(t/60,[m_d1,m_Fd_1],'LineWidth',2); ylabel('Concentration [mM]');
% yyaxis right; plot(t/60,-gradient(m_in,t),':','LineWidth',1.5); ylabel('-dm_{in}/dt [mM/s]');
% xlabel('Time [min]'); xlim([0,min(150,t(end)/60)]); grid on;
% legend('m_{d1} (slow rail, output)','m_{P1} (fast rail)','input falling slope','Location','best');
% set(gca,'FontSize',fontsize);
% 
% % 4. Differentiator 2: the two rails, and the slope of m_d1.
% figure;
% yyaxis left;  plot(t/60,[m_d2,m_N2],'LineWidth',2); ylabel('Concentration [mM]');
% yyaxis right; plot(t/60,gradient(m_d1,t),':','LineWidth',1.5); ylabel('dm_{d1}/dt [mM/s]');
% xlabel('Time [min]'); xlim([0,min(150,t(end)/60)]); grid on;
% legend('m_{d2} (fast rail, output)','m_{N2} (slow rail)','slope of m_{d1}','Location','best');
% set(gca,'FontSize',fontsize);
% 
% % 5. Check of the identity: difference of the rails against the ideal subtraction.
% figure;
% plot(t/60,diff_1,'LineWidth',2); hold on; plot(t/60,ideal_1,'--','LineWidth',1.5); grid on;
% xlabel('Time [min]'); ylabel('[mM]'); xlim([0,min(150,t(end)/60)]);
% legend('m_{P1} - m_{d1}','w_1 (fast-delayed - slow-delayed m_{in})','Location','best'); set(gca,'FontSize',fontsize);
% 
% % 6. Counter.
% figure;
% yyaxis left;  plot(t/60,m_c,'LineWidth',2); ylabel('m_c [mM]');
% yyaxis right; plot(t/60,m_d2,'--','LineWidth',1.5); ylabel('m_{d2} [mM]');
% xlabel('Time [min]'); xlim([0,min(150,t(end)/60)]); grid on;
% legend('m_c (counter)','m_{d2} (consumed by the counter)','Location','best'); set(gca,'FontSize',fontsize);
% 
% % 7. Neuron: output against w_N*m_c + b_N.
% figure;
% plot(t/60,m_out,'LineWidth',2); hold on; plot(t/60,target_N,'--','LineWidth',1.5);
% plot(t/60,b_N*ones(N,1),':','LineWidth',1.2); grid on;
% xlabel('Time [min]'); ylabel('[mM]'); xlim([0,min(150,t(end)/60)]);
% legend('m_{out} (neuron)','w_N m_c + b_N','bias b_N','Location','best'); set(gca,'FontSize',fontsize);
%%
% 8. Enzyme Concentrations Bar Plot
figure;
% Combine the calculated enzyme values in uM

% Create a cell array matching the labels for all enzymes
enzyme_labels = {
    'r_F'; 'r_Fr'; 'r_S'; 'r_Sr'; 'r_I'; 'r_Fx (new)'; ...
    'r_F_2'; 'r_Fr_2'; 'r_S_2'; 'r_Sr_2'; 'r_I_2'; 'r_Sx_2 (new)'; 'r_c'; 'r_I2'; ...
    'r_Nw (Neuron)'; 'r_Nb (Bias)'; 'r_Nr (Removal)'; 'r_I1 (Counter)'; 'Block 0 (Total)'
};

% Create the bar plot
bar(E_all*10^3, 'FaceColor', [0.2 0.6 0.8]);
grid on;
set(gca, 'XTick', 1:length(enzyme_labels), 'XTickLabel', enzyme_labels, 'TickLabelInterpreter', 'none');
xtickangle(45); % Rotate labels so they don't overlap

% Formatting labels and titles
xlabel('Enzyme / Reaction Link');
ylabel('Enzyme Concentration [\muM]');
title('Calculated System Enzyme Concentrations (\muM)');

% Dynamically adjust font size if defined in parameters
if exist('fontsize', 'var')
    set(gca, 'FontSize', fontsize);
end
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
%   data: {"layout":"onright","rightPanelPercent":35.9}
%---
%[output:4c4d511c]
%   data: {"dataType":"text","outputData":{"text":"Adding matlab path to: D:\\code\\bioelectrodynamics\\code\n","truncated":false}}
%---
%[output:96547eb7]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAeAAAAEhCAYAAACnYHOEAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQucjdX6x39zM8aYYWYMMy41JJfklKQI5RSHEkUl5BK6E0pyK7mGipIQJaTcIsWhlEtJhzgqf3VEQmXMMMwwM25z2\/\/Ps7R37+zZe7\/vu\/f77v1envX5+JzTzFrPWs\/vWXu+e93DHA6HA5xYAVaAFWAFWAFWIKgKhDGAg6o3V8YKsAKsACvACggFGMDcEVgBVoAVYAVYgRAowAAOgehcJSvACrACrAArwADmPsAKsAKsACvACoRAAQZwCETnKlkBVoAVYAVYAQYw9wFWgBVgBVgBViAECjCAQyA6V8kKsAKsACvACjCAuQ+wAqwAK8AKsAIhUIABHALRuUpWgBVgBVgBVoABzH2AFWAFWAFWgBUIgQIM4BCIzlWyAqwAK8AKsAIMYO4DrAArwAqwAqxACBRgAIdAdK6SFWAFWAFWgBVgAHMfYAVYAVaAFWAFQqAAAzgEonOVrAArwAqwAqwAA5j7ACvACrACrAArEAIFGMAhEJ2rZAVYAVaAFWAFGMDcB1gBVoAVYAVYgRAoYAsAFxQUYOzYsVi5cqUqiTt16oSpU6ciJibGZ7nMzEysWrUKGzduxIEDB1BUVISkpCQ0bdoU3bp1Q8uWLREdHa2qbj1sqmoAZ2YFWAFWgBXQVQFbADgnJwePPfYY9uzZo0pMOQATaFesWIEJEyagsLDQq+1mzZph4sSJqFevnmz9etiUq\/TYsWOgf5xYAVaAFXBXoGbNmqB\/nLRXwBYAPnjwIPr164eMjAxVCvoCMIFy1qxZ4p8zRUZG4qqrrkK5cuVw9OhR5OXluX5Xu3ZtzJ071yeE9bAp5zCBt+CZWiJb2z215bLz71kBVsBGClQrV4SxHRuh3ci5iExOs5HnwXHVFgDevn07+vTpIxS9+eabMWXKFMTHx8sqTECNi4tDWFhYmbxbtmzBE088IaabKfXo0QNDhw5FcnKy+G\/6+bZt28TU9\/Hjx8XPmjdvjpkzZ7ryuBvVw6ackzt37kTV6S1Etj39t6JWQnm5Ivx7Cyvw3XffiT766quvokaNGhb2lF1TosDF\/32NKz4bidRxWxHTqI2SIpxHhQK2APB7772HSZMmCVl69+6NF198EQRXfxNNaQ8cOBAEL0oE4meffdajTRp9P\/nkkzhy5IjIO27cONeXAWn9ethU4p8UwJHzLuCKRAawEt2smodmbhYvXoy+ffsiLY1HPFaNs1K\/Lvz8FTLG\/ZMBrFQwlfksD2D3DVjTp09Hly5dVMpUOvv69evx9NNPix\/WrVsX8+fP9\/nHitaJR40aJfLTKHj27NlISEgoZVQPm0qclAL4+JQctKpbWUkxzmNRBS5evCiWalJTU1G+PH8Zs2iYFbsVFAA7HMhePQlnVo4DHCWo0OJBpDy7XHEbzZzR8gDOzc3FoEGDQNPQsbGxeP\/999GkSRO\/Y+YO9AEDBmDEiBE+R9T0B43y\/fLLLyIftYFA7Ex62FTqoBTA3z2fgR7NUpQW5XwWVIABbMGgBuBSMAB88cB\/kDGlMxznTouWMoADCJjRitLU78MPP4w\/\/\/wTDRo0wIIFC8S3e39TVlYWHnnkEezbt0+YePfdd3H77bf7NOcOWJqupi8FzqSHTaX+MYCVKmWPfAxge8RZqZd6A7joTCaOv3Q7io7vdzWJAaw0OibIR4Dp2bOnaGnHjh3FuV4aCfubfvrpJ7Gj+vTp0wLkCxcuVHS86K233sKMGTNEtTQF\/vLLL7vOButhU6l\/UgB\/9MhhjGjPO6GVamfFfAxgK0bVf590BXBRAU68\/SjOff1+qQYygP2Pl+FKLl++HKNHjxbtos1S99xzDzZv3oyvv\/4aP\/74ozi\/GxUVheuvvx4dOnTA3Xff7XWXMtmQ7qhWM6KWrgPTBR20buxcB9bDptJASAG87J6tGNOLdzoq1c6K+RjAVoyq\/z7pCeC8r5bg5NsDEFZc+g4FBrD\/8TJUSToKRBdgLFmyRHG7CMbDhg0Tu0A93V4lBel1110nprQTExNl7UtH4u7l9LAp26C\/MjCAlSplj3wMYHvEWamXegH44qHdyJh0p1j3jaxaG+EVE1Fw+PJFSQxgpdExeD7pBiz3ptIZxypVqoDWZw8dOlTmJisCMIG4YsWKpYpKYalmSlsK4Fq1amHRokWgyzko6WFTaWikAJ7WegPmDb5TaVHOZ0EFGMAWDGoALukB4OJzZ3Bi6j24+Ms2OCKiUW3gQuTv\/hTnd6xgAAcQK8MVpY1XtF57+PBh0TYa3T711FNiTdh5YQb9nEbKu3btEpcP7N271+XH8OHDxbS19CIO6ZliuasqpYJIAVy9enWxdnz11VeLLHrYVBoMKYAHNvoYn40L7IiW0no5nzEVYAAbMy6hapXmAHY4cGrJ88hd95pwKa7dE0h+dA4yX+\/BAA5VkPWql3ZA05nb\/fv3ixHutGnTfB5BOnnyJAYPHixgTIkgTSPVhg0buppINubNmyf+WysA62FTqaZSAHdKXYo9rzCAlWpnxXwEYHoIpFq1arKPkFjRf\/aptAKFB7Yja1I7zS7iyP9mGU7M7oew4ksoV68FUkdvQERsZWTO6M4A5s4H0I5kOrNLR4MoPfPMM+LIkHMULN3UpRWA9bCpNJbuAP73M3+fT1Zqg\/NZR4FLly6Jvk9fPtW+4GUdFdgTlwK\/7YRjXk9NAFyU8RuOT2qHopNHEBabhNRRa1G+\/i2iKgYw9zmhgPuZ3bZt24rjQ861YOl6rb8A9rUGrJVNpeF0B\/CakXfydZRKxbNgvgsXLuDEiRNISUnhm7AsGF+1LtEUdM7UOwMGcMn5XGRMuweX\/vcVEBaOxO4TUbnr5dMpDGC1UbF4fumItHHjxuKyDeeasR47lvWwqTRE7gCmTVh8HaVS9ayXj9eArRfTQDzSZA3Y7arJ2JY9UG3QIiCyHAPY4XA4AgmQFctKgei+YUq6mcodzr60kNps1aoV6GIO54tMethUGhd3AL\/Q+598HaVS8SyYjwFswaAG4JIWAD7\/f1tw4pXOcFw6h4iUuqgxZiMiU+qUahVPQQcQJCMXzc\/PF9PK9Eav+5Eib+2WwtL9sg3p28LucPalg\/QmrG7dumHChAmiTZT0sKk0Ju4A7tOxOd+GpVQ8C+ZjAFswqAG4FCiAizIPI31yexRnHhJHjlKGfIDYFveXaREDOIAgGbEorWWNHDkS69atE81T8wyhFJbuo1X3e5vpYQXK4yvRxha6jWvNmjUi2wsvvID+\/fu7iuhhU2lM3AHcpsX1mN3j713fSu1wPmsowAC2Rhy18KIo6ygIwFmz+\/m3Bux21WR8p+dQpfcrgIf31RnAWkTMQDZoZn3KlCli\/ZYSXTVJ1z\/S5Ru+Eo2Y6bGETZs2iWzurx25b9Jy3yXtybba15C0sKk0FAxgpUrZIx8D2B5xVuJlzsrxyPlonMiaOm4rYhqpu6b27PqZOL34WfHEYPkGt6LayE\/FkSNPiQGsJCImy0P3PT\/66KPiog16BvCNN97AXXfd5dOLLVu2iMs3nGXeeecd3HbbbaXKSN\/upTPCBHa6WctbUvsesFY2lYTLHcBRVdPw4wstlBTlPBZUgAFswaD66RKNgPO2LhYQVgtg6ROD4sjRC5+hfN1mXlvCAPYzSEYulpOTg4EDB4IgQ4mufpw7d67X14toLfbJJ58EXeBByduRIJoyprPBu3fvFvnodaPx48d7XGN2tzlu3Dj06dOnjGx62FQSG3cAZ0SmIHvGP5UU5TwWVIABbMGgBuCSv2vAOR9NQM7KlwKo+e+iVYetQsXm92liy2hGwqy+C5pGwQTLc+fOCe2vvPJKjB07Voxqw8PDxc9ovXjt2rV47bXXxDODSmC9YcMGDB06VIyUnRCmqyvp\/CSlkpISfPXVV6Ku48ePi581b94cM2fO9Prakh425TqcFMCPV52BPdHXixHwFYnl5Yry7y2oAAPYgkENwCUGcADiKShqeQATCD\/66CO8+OKLLliSLklJSahZs6bHxxjod6+\/\/rrPzVUE3lmzZol\/0kSjbDpedOzYMRfMlQCd8uhhU64PeAIwbcLq0ezyFwlO9lKAAWyveMt5ywCWUyiw31sewCQPbciiN4BHjRpVCoqepGvWrJl4wrBevXqyytLu5sWLF2P69OllXlOSFg61TV+OMIBlw2yrDAxgW4Vb1ll\/ASxr2C0DrwGrVcyE+Wka+ssvv8Tnn3+OPXv2lJpupqNEtEGLYOmcmlbqIo12V61aJSB\/4MABMZKlUXTTpk1BZ35btmyp+l5dPWx68scTgGn0y0eRlEbfWvkYwNaKZ6DeMIADVdB3eVuMgPWV0NzWPQG45VWVsW5gE3M7xq33SwEGsF+yWbYQA1jf0DKA9dXX8NY9AZgazTuhDR86XRrIANZFVtMaZQDrGzoGsL76Gt66NwDzRizDh06XBjKAdZHVtEYZwPqGjgGsr76Gt+4NwLwObPjQ6dJABrAusprWqBUBXJJ7Crmb5iP3mw9RmPEbwooviXuqo1KvQmzTu1H5rsGISPR+sZKWwWQAa6mmCW15AzCdA+YbsUwY0ACbzAAOUECLFbcUgIsKkP3JKzj78RQ4Cs97j1R4BOLaP4UqD01FWHQFXSPKANZVXuMb9wZgajlfyGH8+GndQgaw1oqa216wAKy7Sm4PQyipT+7+aiU25PIwgOUUsvjvpQCeXH8u1pyv7\/KYp6EtHnwP7jGA7RdzXx5bBcDShyHURDj2tj6o9sQ7QOTlp2O1TgxgrRU1mT0pgL\/quAzP\/V+1Uh7wbmiTBTTA5jKAAxTQYsWtAOCinAxkjG2DwsyDruiEV0xEfMehiG\/1EMLKRcNRUowL\/\/sGZz+dhoI\/9rnyhUXHotrza1HhH7frElkGsC6ymseoFMD04kmLNdH4I\/uiy4ER7dMwon1t8zjELQ1IAQZwQPJZrrAVAJy7dTFOze0vnkWkVK5eC6Q+uwoRSdXLxMtx6Tyy5j2O\/G8+cP2u4h2PoOrj8z2+YxxowBnAgSpo8vLuAH7z2JWYtvFoKa\/WPtUErep6fsfT5O5z890UYABzl5AqYAUAZ739OPI2zxduhcdVReqLXyC69nVeA027pNMntkfh0e9FnqjUBkid8DUiK1fVvHMwgDWX1FwGpQBOHrgQcW0exvWTdpQaBfOOaHPFNJDWMoADUc96Zc0OYMfFfGRM7oiLv2wTwalwQyekDF8NREb5DNbZtdNxeslzl6Edm4DUsV8iuk5TzQPMANZcUnMZ9ATg7YfOoPOcH0o5whA2V1z9bS0D2F\/lrFnO7AAuOX8Wxyd2QMGhy2\/CJ3Qbj4QHxsoG68K+rcic2hGOggtwREQhZcSniG1yp2w5tRkYwGoVs1h+TwAmFzf8dAq93vt7MwL9jCBM09H8VrDFOoHEHQawdWPrj2dmB3BxdjrSx9yColN\/qAJwYcavOP7ibSg+myHKVR22ChWb3+ePhD7LMIA1l9RcBr0BmLy4Z+6P+ObXnDIO8cYsc8VYTWsZwGrUsn5eqwE4ecgyxLXqLhs4d3AzgGUl4wz+KOALwGTvgfl7sfmX7DKmOzSqgqldrubRsD+iG7gMA9jAwQlB06wGYKUgZQCHoLPZsUo5AJMm0zYeKbMz2qkVj4at1WsYwNaKZ6DeMIAvT10rBbdavXkKWq1iFsuvBMDksqeNWU4paE34re4N+aiSBfoGA9gCQdTQBQYwA1jD7sSm3BVQCmAqRxd0DFy2H9\/+dsajkDwaNn\/\/YgCbP4ZaesAAZgBr2Z\/YlpsCagDsLLpsd6YAsadEo2G6OYvukeZkPgUYwOaLmZ4tZgAzgPXsX7a37Q+AnaPhZbszvK4Nt2+UhGld6vEmLZP1MAawyQKmc3MZwAxgnbuYvc37C2CnarQ2PGj5\/lI3Z0kV5Wlpc\/UvBrC54qV3axnADGC9+5it7QcKYKd4vnZK8yYt83QxBrB5YhWMljKAGcDB6Ge2rUMrADunpX1t0nq+fRpG8stKhu5rDGBDhyfojbMagIMlYGSVK1Bj8n8QkVjDZ5V8DClYETFoPVoC2Omi3CYt2qDFTxwas0MwgI0Zl1C1igHsn\/IMYP90s10pKYATHhiHhG4vaaIBHVnytUnrquQYrH78et6kpYna2hlhAGunpRUsMYD9iyID2D\/dbFdKLwA7hSQQ08tK9L+eEm\/SMlaXYwAbKx6hbg0D2L8IMID90812pfQGsFNQ3qRljq7FADZHnILVSrMDOFg6+VsPrwH7q5xFygULwCSX3E1a3ZpWw9sPXWMRZc3pBgPYnHHTq9UMYL2UvWyXAayvvoa3HkwAO8WQu1eaN2mFrtswgEOnvRFrZgDrGxWfAM7Pzwf9gb755psRFxenb0vYekgUCAWAnaNhX5u00pJi8MmTvEkr2J2CARxsxY1dHwNY3\/j4BHB2djYGDBiA\/\/3vf+jSpQv69u2L+vXrIzw8XN9WsfWgKRAqADsdlJuW5k1aQesKoiIGcHD1NnptDGB9I6QIwHv37nW1ggBMIL7jjjuQnJysb+vYuu4KhBrATgd5k5buoVZUAQNYkUy2ycQA1jfUqgHsbE5kZCTatGmD\/v3748YbbwT9NyfzKWAUADunpX3dpHVfk2p4pzdv0tKzlzGA9VTXfLbNDuCiM5nInNYFRdnpQRU\/MrEGUkasQWRl36\/C+QTwpUuX8OGHH2LBggXIyMjw6kBqaqoYFdM0NY+KgxrngCszEoCdzvAmrYDD6rcBBrDf0lmyoNkBXJydjvQxt6Do1OU7nYOVND0HXFJSggMHDmDx4sX47LPPkJeX59EPHhUHK7za1WNEAPO0tHbxVWuJAaxWMWvnZwD7F19NASxtwoULF\/DNN99g4cKF2LNnD4qKijy2kHZN9+zZEz169ECtWrUQFhbmnydcSlcFjAxgclxuk9bsHg1Bx5Y4aaMAA1gbHa1ihQHsXyR1A7C0OVlZWfj3v\/+NlStXihGyt9SsWTOxm7p169aIiYnxzyMupYsCRgcwj4Z1CbtXowzg4Opt9NosCeCwcJSr1QiV7hyC2Fu6IbxC6I7YanIRh3OK+qOPPsLq1au9TlHzqNh4HzezANg5Gqbd0vTakqfU6R\/JmNi5Lj\/wEEA3YwAHIJ4Fi5odwCX5OTgx+2Fc+H49UFJcNkLhEYhu0AqVOwxEhaZ3I6xccAeImgDYfYp6x44dYvMWTVV7m6KmUfGDDz6Idu3a8SUfIfzgSgEc1+ZhJA9cGMLWKKuaNmkNWr6fH3hQJpeqXAxgVXJZPrPZAewMUMn5PJz771rkrp+JS0e\/9whjR0Q0KjS6FZW7jkZM\/VuAyHK6x1dzALtPUW\/evFls3vI2RR0VFcWXfOgeZu8VmBHAPC2tX4dhAOunrRktWwXAUu1Lck8h95ulyP3sLRSd\/A1wlJQJTVhUBVS46R7Etx+ImPrNgfAIXcKnK4Bd3z5KSvDrr7\/ik08+wdKlS71OUdMlH3SuuEOHDjwq1iXcZY2aGcDOaWlfZ4d5k5a6jsQAVqeX1XNbEcDSmBXlZCB\/60Kc3bQAxVmHPYczJh5xrXqi8l1PI6pGQ0DDDcVBAbDUq4KCAnz\/\/fd477338NVXX3mcor7uuuvE2ePExESr9++Q+2d2ACsdDfMDD8q6GgNYmU52yWV1AEvjWHjyKPK2LETu1gUo8XJxR0Slaoi9tRcqt3sSkSl1AoZx0AEsdZjOE3\/77bf4+OOP8fXXX6OwsFD8mgEcvI+3VQDsHA37euChbcMkvHZfPd6k5aN7MYCD99kzQ012ArArHg4HijIP48yXc3Fu2wcoPnvCY6gikuugUtsBqPjPfohMSPUrnCEFsMPhwPHjx\/H555+L6emff\/6ZAexXGP0vZCUAO1XgTVr+9wcGsP\/aWbGkLQEsDWRJMQqO\/h\/ObHgT+bs+Bi7klg1zWDgiq16F+DsHIb51T4THV1HcFYIOYILuiRMnsHbtWqxYsQJHjhwp01geASuOX8AZrQhgpdPSb3VviFZ1KwesoZUMMICtFM3AfbE9gN1gfOHATuRunI3zuz6Fo\/C8RxhHXXk9Kt8zHHHNu8rupA4agGm6+csvvxTQ3b17t8ee4bzKki7toGNK\/Oxh4B8gOQtWBrBzWtrXJq3n2qVh9J215WSyze8ZwLYJtSJHGcBeZCoqwMVDu3D28zkeYRyUm7DkIkjXVv73v\/\/FkiVLvG64Ihv8xKGckvr93uoAlk5Ld57zg0chr0gsL66zHNGeQcwA1u+zZkbLDGCZqJUU48LP23D6g+dRcPi\/rswhAzBdvLFv3z4x0l2zZo1rY5W7G3z+1xgfR7sA2Dka9rVJi0Bs92lpBrAxPpdGaQUD2EMkZEa\/VCKoAKarKH\/\/\/XfQVZT07\/Tp0177D9+AZZSP1uV22AnATuXlHngY0T7NtqNhBrCxPp+hbg0D+K8IKICuM1bhiTUQf1sfJNz3AsKiK\/gMYUBrwEofY+A7oEP9MfJevx0B7FSD7pWetvGo12lpmpK220tLDGDjflZD0TJbA1gFdP09H6wawLSZatu2bXj\/\/fe9bqYSQ\/DISDRt2hT9+vXjV5BC8clRWKedAeyclva1SatDoyqY2uVq25wdZgAr\/ODYJJvtAExrunI7nf+KfVhsEire8kBAN2QpArDz9iqCLt3t7Lwww1MfTE1NRd++fcX9zsnJyTbppuZ10+4AdkaOzw5fVoIBbN7Psh4ttwWAVUAXMfGoeFNXVL5rMMql\/SPgO6J9Ajg\/Px9vv\/222MVMI19viTZU3XHHHejTpw9uvPFGMfrlZA4FpAAu36gNqo\/bao6G69RKuWlpq2\/SYgDr1LFMalZTADscuPTrLpz94m1c\/HkLiui6R3oiMCwcEfHJiK5\/Cyp3HILyDVoFDDZZuf2Abtwdj2j+MINPAGdnZ4PO5O7du9ejP\/x4gmyYDZ+BAVw2RHKbtKx8dpgBbPiPbFAbqBWAizJ+w4k3e+LSoV2y7Y9KuwEpQ5YgquY1snlVZfgLunmb3\/V+q5VzetkIryF5AjAfH1IVcsNnZgB7D9Gy3Zmg9WFPyapnhxnAhv\/IBrWBWgD44sGdyJx6D0ryTpZqe1hcFYRFxcBReAGOvFOlfxebhNRRa1Ge3uUNJKmArngPuMmdqHz3EGO8BywFMB0fotFw69atERMTE4gkXNZACjCAfQeDRsN2OjvMADbQh9MATQkUwPT2bvrE9ig8+v1lb8IjEH\/n00jsOqbUncmU79TS0cjf+t7laWkAUSn1kPrSl+JMrb+pODsd6WNuQdGpPzyaENBtdCsq3TUEFa5rJ3t1pL\/t8FbO5xR0Tk6OuEyjbdu2qFWrFsI0fAdRa0fYnn8KMICV6SY3LW2Vd4cZwMr6g11yBQrgs+tn4vTiZ8Wj9wS7qk+8g7jbenl+xs\/hwJl\/v4HsJc+J\/JQSuo1HwgNj\/ZbbE4Bd0L3zacT8oy3CyoVuQKloF7Tf3nNBwyvAAFYXIrlNWma\/0pIBrK4\/WD13IAB2XDqPjKmdcPGnLUKmmBvvQeqwlT5Hme5louu1ROqYDQivEO+X1HIjYL+MKigU1JuwFLSHsxhUAQaw+sDQaJhATGvEntLt9RMx44H6pjw7zABW3x+sXCIQABdm\/IrjL96G4rMZQqIqTyxA\/B39ZeU6s\/5NZC8aIvIpBZk3owxgWbkvZ\/jhhx8wYcIE8f9r1qwp\/n9CQoLC0pzNXwUYwP4qB1jx7DAD2P\/+YMWSgQC44PhBZM3qjcLTx4GCc6g2bDViGv9TVqb8natxcvr9Il94bAJSx36J6DpNZct5ymBqALvvgl66dCmaN2\/ulxByhQgEPXv2FNn4PWA5tbT7PQM4cC3lpqXNdHaYARx4f7CShUAA7K8OOR9PQc6y0ZqMgIvOZCJzWpfLZ46DmCITayBlxBpEVk7xWauqc8BKAOwvtBnAQewdkqoYwNroLrdJ6+7GyZh0T13DT0szgLXpD1axEmwAa70GbPQ4MICNHiGd28cA1lZgmpb29u4w1WT0l5YYwNr2B7NbCzaAz+1YhcyZvRBWfElIl9T7NVTqPMzsMnptPwPYsqFV5hgDWJlOanOZdVqaAaw20tbOH0wAF2UeRvrk9ijOPCREjah6Faq\/tAlRVdNMJ7Lj4jmElSsve6UmA9h0odW2wQxgbfWUWpOblr7\/hmqY30vj6\/YCdIcBHKCAFiseLAAXnzuDE1PvwcVftgkF6axutYELUbF1D\/Mo6nCg4Pd9yP54Mgp+240aE79BRGINn+1nAJsnvLq0lAGsi6yljMqNho10dpgBrH9\/MEsNRVlHQQDOmt0PqeO2IqZRG12aLuD72v24+NPmy\/bDwlHx9v6o+sgcIDJKlzq1NFpyPg9525fizCevoDjrsDCt9PgUA1jLSJjQlhTAkclpuGLOERN6Yfwmm+VKSwaw8ftSsFqYs3I8cj4aJ6rTC8DFp48jY8b9KDi4w+VW7G19UO2Jd4J+LaRaXS\/98TNyVk3EuV2fuNasnTYYwGrVtGl+BnBwAy83Lf38v9IwskPt4DZKUhsDOGTSG65iGgHnbV0sIKwHgAuOH7h8ROj4Xw+e0Mj31j5IfnQ2wqIrGE4P0aCiAuT\/dx1yVr+Mwt9\/dF2Z6d5YBrAxw2e4VjGAQxMSo05LM4BD0x+MWqtea8AXftyIE2\/2+fuFpLBwxP3rKSQ\/PN2QI9\/Ck0dxdu1ryN22BLiQ6zNcEcl1kNh1FOLa9JH1haegjdrzg9QuBnCQhPZQjRGnpRnAoesPRqxZcwA7HDi7bgZOLR3jmralDVdVek5GpU7Pen6kIVTClBTj4i\/bcXrxcFw6ssfraJeaRz7E3nQvEruORrkrGyv2gwEcquAapF4GcOgDQSCms8P0v55SMF9aYgCHvj8YqQWaArioAFmLhiHvizkumIVXqITkJ95FbPP7FENLb33oacQzX8xF3uezUXz2hOxot\/K9zyOuVU+EV4hT3TQGsGrJrFWO6ba9AAAgAElEQVSAAWyceBphWpoBbJz+YISWaAVguuEqa97jyP\/mA5dbkVVrI2XEpyh3RePQu\/rXaDf7o4m4+L+vXG8Se2xYWDgqNLsXiQ+8pGq068kWAzj0oQ9pCxjAIZW\/TOVyLy1dkVgeet4tzQA2Vn8IdWs0AXBRAU68\/SjOff2+y51y9Vog9dlViEiqHlIXPR0hkmsQjdpTx2xEdL2b5bLK\/p4BLCuRtTMwgI0ZX7mXlvSalmYAG7M\/hKpVAQPY4UD26kk4s3Kca9o5tsUDSH5igV9Ttpro4Low42WPR4ikdURUqobYW3shtlkXZL3ZE0Wn\/gADWJMosBFSgAFs7H4Q7GlpBrCx+0OwWxcogN3vdi7f4FZUG\/kpImIrB9sVcYQob\/tyZK8c77oww1MjnBuqKt35NGLqNxfXSUqfNWQABz90lq2RAWz80MqdHb4yqTxmPdgQreoG\/keNAWz8\/hDMFgYCYNrMlD6xPQqPfi+aHE5P9A1dhqhqdRS74AgLR2RckuxxHl8GFR0hCgtH1JXXo3KnZxF7Y+cyo3MGsOKQqc948uRJDB48GLt27UL16tWxcOFCXH311bKG1q9fj6efflo2nzSDUvuZmZlYtWoVNm7ciAMHDqCoqAhJSUlo2rQpunXrhpYtWyI6OlpV3Z4yM4ADljBoBoIxLc0ADlo4TVFRIADO3boYp+b293l8R06EQEebUnB6qou+FMS3exzxdzyCyIRUr80xBIDnz5+PG264wadmZ8+eFTD7+eefRT4lZSjf999\/j8cee0yUue6667BgwQIkJibKxSfg3xPYZs2aJf5RUgpIyjt9+nTMnj1bVRvk7FN7VqxYgQkTJqCwsNCr7WbNmmHixImoV6+eqvrdMzOAA5IvJIX1nJZmAIckpIat1G8AOxw4MaMbzu1cFZBvegCYoFuxZXdUbvckIlPqKDr+ZAgAB6SkisLBBPCGDRswdOhQMcJUA+Bz585h5MiRoFGwmuQLwO5fBshuZGQkrrrqKpQrVw5Hjx5FXl6eq7ratWtj7ty5AUGYAawmesbJq9e0NAPYODE2Qkv8BXDJ+VxkTL4Llw5+G5AbWgM4uu5NiG31ECre0FExfMkBBnBAYfRc+ODBg3jyySdx5MjfDxDIjVCdlrKysvDII49g3759iI2NxZw5c9CoUSPZVoaFhSEuLk6A1T1t2bIFTzzxhOvLQI8ePcSXg+TkZJGVAL1t2zaMHTsWx48fFz9r3rw5Zs6c6coj2wC3DFIA06\/qfORQa4Lzh1ABraelGcAhDKYBq\/YXwEZxxdcUdFhUBVS46R7Etx\/o2mzlrd0MYI0jmp+fj1GjRpUZwSoF8E8\/\/YR+\/frh9OnTaNCggZgyT031voYg1\/ycnBwMHDhQ7EqmRCB+9tlnPYLa\/YvDuHHj0KdPH7kqPP6eAeyXbIYrpNW0NAPYcKENaYPMDmC6AOTs+jdwdsObvm+1iolH\/K29UbnjMx5HxiEBMEFq5cqVoHXdYKZKlSqJjUYVK1bUpVoaSdLaLY0cKXXv3h1ff\/01MjIyFK8BSzdgdenSBS+\/\/HJAm6Kk9urWrSvWztPS0rz6T+vE9AXCOQomfxISElTrxQBWLZlhC8hNSyu5xIMBbNjwhqRhZgewS7S\/brrKWfcGzv\/wWZnnA6Xi0mMKlToOQXzrngiPryJ+FRIAhyTiQah0+\/btYrRJ66kEz969e2PQoEFiWlfpCPitt97CjBkzRGtppErl\/U0FBQViWpm+7FAaMGAARowY4XH066yDvixQvl9++UXke\/\/998V0tNrEAFarmPHzBzItzQA2fnyD2ULLAFgimuLbr8IjEJ12Ayp1GYHoOjcg46U2wb2II5iBDlZd6enpYrf1\/v374Rxp0m5jmk5WCuBLly5h9OjRWLNmjWj2u+++i9tvv91vF6TryUrtuUPb3y8BDGC\/w2b4gv5MSzOADR\/WoDbQigB2CehwoCjzMM6sf132mUFHWCTCUCKOVAW6MUwaQJ9XUQY10kGojKbUX3rpJQFO2jhFo9jbbrsNv\/76qyoAZ2dni9Hn3r17xbovnRsO5DiQdD1ZjT3pKNzfaXAGcBA6XgirUDstzQAOYbAMWLWlASzVu6gA+f9dh5zVL6Pw9x9lnh6MQlybh5HQdTSiqnpfJlQSTtsA2OFwiCle57opXaBB\/2j6Vi2AaRMUjZhpGpguxpg8eTIIol988QX27NkjNmZRomNCbdq0wf3334\/69esjPDzcY0xoSty5iUrNhi7pOjC1g9aN1a4DM4CVfEzMn0duWvqt7g3Q86ZUMIDNH2stPbANgKVT1EqfI6RbuqpehUqdh+nzHKGWgQy1LQLN448\/LtZ93Y\/uqAUwHReiI0hq0l133YUXXngBKSkpZYpJQarmDDT51LNnT2FPTTlpAxjAaqJo\/rxy09IPNKmCntdEipmd8uXLm99h9iAgBewIYJdgKjZu0f3RFRrdioT7xqB8g1bi\/mglyRYjYPerJunIEI1InUktgKVTv1KR6Xwv7Vymke6xY8dcI2FnHoIk3Z5Vp07pu1ClAO7YsSOmTp0qpsjlkhTAtWrVwqJFi8SoW01iAKtRyxp55Z48rB4fKZ48bNPg8g5QTvZVwNYAloRd8cYtAGGxSYhv0xeJPSYiLLqCz85jeQDTkSM6IkRwoulmmi6mKWG6EMMfALtvfiIbbdu2xTPPPFNqmpmmvP\/8809xQYdzdzPlbdWqlTj+JJ0qfu+99zBp0iTRnE6dOgkAx8TEyH7qpQBWunvb3ag7gK+YcwSRyYGta8g2nDMYQgG5aelh7a7EmDuVX5xvCKe4EZoqwAB2k9PhQGH6fpzZMAt525cCF3I96h1Z5QrUmPwfRCTWsC+ACYL0oMGYMWPELVJdu3YV9ye7w03NCJiuoKRNV3Ru+I8\/\/hDHmWga2NPNVqS8p+sl3S\/OmDZtGubNm2cIACdP3x\/wxgJN\/wKwMd0VoGnpmV9neqyHzg73aJaCIbf5f8mM7g5wBbopUHhgO7ImtUPquK2IadRGt3pMadjHxi0GMADpjVENGzYUm5Rq1Cj7jUQNgP3pKLTuPGzYMGzatEkUpzVo6cUZy5cvF8eajDACDhu1DUio6Y+bXMbEChzPLcK6\/fmYv+uMRy9oWnpc2ypoWoPXhU0cZvVN\/20nHPN6MoBllKOnF898MRd5n88WN27ZHsB0teOQIUNAO4xpbZaAR9O\/npLeAKY6pTdd0bOCNIq+9tprRXOka8D+TkFrtQac\/MKXiKrvWSf1n14uYTYFDp\/MQ7f3fgEB2VO665oEjL0rDTQy5mR9BWgKOmfqnQxgpaH+a+NW7pZFSOw1BZGVy266lZqy5Bqw+7QvgZimir1NEwcDwNKjSxQAurnK+YXASLugeapJ6SfNmvmcx5CW\/q8I0zcf8+rkiPZpGNFe3YY\/aypmba9suQb81wUdjuICRNVs6DXAZz6ZhsKTRxDfph+ir75J0bOG7sYsCWDpRRmBfjzUjEh91SWFPOWbMmUKHnzwQVFEupmqcePG4mYt5wtIvmxKwU0wp93Z8fHxqlx234TFAFYln+UyS88BnzwPLNudgWkbj3r0U8nd0pYTyGYO2QrADgcu7P0CpxYMRmHmQcTd8RiSn7i8N8dTypzRHed3rBC\/ikq7AdUGLkS5tH+o6iEMYBm5PAGYdkLTrVq0k5oejvB2wYbUtDuApddXSkfHanYzS49D0eMVEyZMEO8Gq0nuAE4euFDc8sLJngp4uohD7jatOxokYvr99Xla2oJdxjYAdjhwZs0UZK8YC5QUi0jG\/KMdUkZ8irByZU+keHrvmJ43rPLUAsS16q64JzCAVQJYOupUc2vVDz\/8IG67ol3UdMaXpqCbNGkiane\/C1o6Pe2tee73UdMlH\/3791cceGdGBrBqySxdwNdNWL4u8SBRZvdoKHZMc7KOAnYBcP43y3Bidr9SryRFpTVBjZc2IbxiYpmAFp05iYzxd6Dw2E+lfkdngFNHrUX5+rco6gSWBLAizyWZ1KwBS6+NVPMKkfSsL00z07Ej561Y7meL6Uwxva4kPavs7pNeryHxCFht77FWfrmrKGk0LDctTRDm9WFr9As7ALgoJwMZY9uIaWeRwsJR\/rp\/IaHrC4ip39znrVaFJ48iZ9UE5H\/9vmvkHF2vJVLHbEB4BfnlQAYwoOou6MzMTHGl5b59+0SsevToIZ4SjI6O9vqJk77ARJk8PTco3SXt68iUsxK93gNOeGAcErq9ZI2\/HuyFagXkAOw0KDctzevDqqU3ZAE7ADh383s49fYAoT9dKVml\/0xUaveY8k1VDgfyvv4AJ99+VIygHRFRSH1uNSrc2Ek2pgxglQCmyz1o7fX1118X4nq7XcupvPQFJvqZtzVemoamUe\/u3btFUXrdaPz48ahYsWKZIErPN9Mv3S\/2kI26JIP7FDQDWI161surFMBOz+Vu03I+8mA9pezhkeUBXFSIzFfvw\/nv14mAVmzdC1Wfeg+IjFIX4KJCnJzTH\/nffCDKxbbsgWpDPpSFOANYJYBJXBoF00tK9PIRpaioKPHf9EKS8w7nkpIS8VwhXSvphKocrDds2IChQ4eK27OcEB4+fLhrqppsfvXVV2LETW8XU3J\/WEJdr7m8A7vq9BauYgxgtQpaK79aADu9l3vkgaelzdlPrA7gkvNncXxiBxQc2ik2W6WMXI+Yxv\/0K1gX9m1F5tSOcBRcQLm6zVH9xc\/F28G+EgPYDwCToPv37xcj1iNHjrj0JRDXrVtX7ER2f4yB4EtnkX2dR\/Z0bSUZpwcW6HiRu036+dy5cwN6i9gdwLQDmtaBOdlTAX8BTGrJPfLA09Lm61NWB3BxdjrSx9yColN\/KL69ylsUycax0S1RknNMsS0GsJ8ApiDQYwsjRowQo0hfiW6+ovuoO3fuLHtkiXY3L168WLyaVFhY6NVss2bNxL3W9erVC+hTzQAOSD7LFQ4EwE4x5Kal726cjEn31OVjSyboPQxg5UHyB+YM4AAATKGhaWGaYqbpY9oh7RwRE3SbNm2KDh06oF27doqeF5SGmka79JDE5s2bceDAATEt7bRJZ35btmzpc+OX0m7DAFaqlD3yaQFgJdPSlIePLRm\/T1kdwOI40djbUJjxC8JjE5A69ktE12nqV2AK\/tiH9PF3wJGbhciqtVF90reITPD9iAkD2C+prVOIAWydWGrhiZYAdk5Lyx1boiNLfH5Yi+hpb8PqAKb12sxp9+DC\/30pxEvq\/RoqdR7ml5C5Wxfj1Nz+gKME5RvcitQx6xFWvuwmWqlxBrBfUlunkDuAyzdqg+rjtlrHQfZElQJaA9hZudyxpSuTymPWgw3Rqm5lVe3lzPoqYHUAk3qnl4zE2bXThJBRKfWQ+tKXYg1XTaLXkNIntkfh0e9FsUodn0XSw9NlTTCAZSWydgYGsLXjq9Y7vQDsbMey3ZkYuGy\/12bxsSW1EdM3vx0AfOnojzg+ri0c504LMaOv+SdSnvsIEXFJisQtOZ+HrLcH4NyOj0T+sOhYpI75DOUbtpYtzwCWlcjaGdwBHJmchivm\/L2z29res3fuCugNYGd9fGzJHH3PDgCG2xleigyNgCs\/MB5xrR70eBe0iF5RAfK\/+wSnPxyF4qzDroCKM8CDFgGR8vfyM4DN8TnQrZUMYN2kNaXhYAGYxJGbluZjS6HvQrYAMLE08zDSJ7dHceah0qKHhSMiPhlRqfURUS1N\/K74xFEUnDiEkrMnXNdPOgtFpNRFjTEbEZlSR1HwGMCKZLJuJgawdWPrj2fBBLCzfXLHljo0qoKpXa7mY0v+BDTAMnYBMMl06dB\/kfFqV5Rk\/+mXauGJtZA6\/GNE171RcXkGsGKprJnRHcDkZZ2PHNZ0lr2SVSAUAHY2il9bkg1P0DPYCcBidHv6OE6++xQufP\/vMqNbr+KHhaNCkzuR\/Nh8RCRVVxUjBrAquayXmQFsvZgG4lEoAUzt5teWAome9mXtBmCngvTKUd62JaBnCovolaS\/3gh2KRwWjsjEWqhwy\/2o1GEQoqpenp5WmxjAahWzWH5PAKZNWLQZi5P9FAg1gJ2K8\/qwMfqeXQHsrn7J+Vw4LuaJH4eVq4DwigmaBIgBrImM5jVCAC6Z0hop5S4\/AEGJ3wQ2bzwDbblRAOz0g9aHO8\/5watbQ++4AmM7XhWo21zeiwIMYH27BgNYX30Nb50A\/MPodmiflO9qK7+IZPiw6dZAowHY6SivD+sWcp+GGcD66s4A1ldfw1snAC96ujOeT8tytZVfRDJ82HRroFEBTA4reW2Jr7XUtmswgLXV090aA1hffQ1vnQA87fGumF4vo1RbeSe04UOnSwONDGCnw7w+rEvoPRplAOurNQNYX30Nb50A3LNnT2xqyrdfGT5Y3EBLKGCm2+YYwPp2OQawvvoa3roTwCvalkNSzgHDt5cbyAqYXQEGsNkjqF37GcDaaWlKS04AL126FFd8PgoXf\/7KlH5wo1kBsyjAADZLpPRvJwNYf40NXYMUwM2bN8eFfVuQNXcAirKOGrrd3DhWwKwKMIDNGjnt280A1l5TU1l0B7CpGs+N1VwBM2zCUuK03LGlNx9sgF43pyoxZes8vAasb\/gZwPrqa3jrDGDDhyioDbQKgEk0uWstKc\/sHg3Ro1lKUDU2U2UMYH2jxQDWV1\/DW2cAGz5EQW2glQDsFE7JsSU+P+y5mzGA9f34MYD11dfw1hnAhg9RUBtoRQA7BZR79pDfHy7b1RjA+n78GMD66mt46wxgw4coqA20MoCdQsqtD\/P90n93OQawvh8\/BrC++hreOgPY8CEKagPtAGBeH1bepbQGcEnuKeSsfRX5\/1mF4tO\/i2f+HBHRiEq9CvGtH0J8h4EIr1BJeQNNnpMBbPIABtp8BnCgClqrvF0ALF0fphHxst2ZHgNJ09K0SYvWiO2YNAOww4HcTe\/i9MKhcBSe9ypleMVEJD\/xDmJv6gKEhVlecgaw5UPs20EGsM07gJv7dgOwFMT07CFt2PKU7Lo+rAmAHQ5kr56EMyvHAY4Sl7xhcVUQFhUDR97pUlCmEXG1gQtRsXUPy384GcCWDzED2OYhVuW+XQHsFEluffhf1yThla71QEC2Q9ICwOf\/bwtOvNIZjkvnhGTlrmiMqk9\/gHJp\/7gsYUkxzv+4ESfnPYGS7D\/Fj8Jik5D6wmcoX7eZpWVmAFs6vPLO8QhYXiM75bA7gJWC2C7nhwMFsOPSeWRM7YSLP225DN96LZA6egMiYiuX+VgVZR5G+uT2KM48JH5X4aauSHlmORAZZdmPIAPYsqFV5hgDWJlOdsnFAP470nIXedAo2OrnhwMF8LkfPkPmtHsQVlwoNlulPv8xKtxwl9ePU9725Tj5Vh+Rn0bB1cdtQnTa9Zb9+DGALRtaZY4xgJXpZJdcDOCykSYQD1q+H3SO2FOy8katQAGc9fbjyNs8\/\/Lot86NSB37pcfRr1PXkvxspI9vi8KjP4gfJXQbj4QHxlr248cAtmxolTnGAFamk11yMYC9R9qOF3kEAuCS82dxfGIHFBzaKUSt1PFZJD08XfajJIV2zD\/aIWXEpwgrFyNbzowZGMBmjJqGbWYAayimBUwxgOWDKLdRq\/01SZhmkY1agQC4MONXHH\/xNhSfzRCiJg9ZhrhW3WUFPrP+TWQvGiLyRVatjeqTvkVkgjUfzmAAy3YHa2dgAFs7vmq9YwArU0xufZiszOreAA\/dZG5wBALgSwd3IH1SB+BCrhjBpoxcj5jG\/5QVOH\/napycfr\/IFx6bIKato+s0lS1nxgwMYDNGTcM2M4A1FNMCphjA6oJIIPZ1kQdZM\/OO6UAAXAqkFSohdcxGRNe7WVbgSwe\/Q8bk9qApbLoVS2k5WcMGzMAANmBQgtkkBnAw1TZ+XQxg\/2IkNmot24\/tv3nfqGXGHdOaATihJmq+\/C0iq1whK7AUwAgLR8rwj1Gh2T2y5cyYgQFsxqhp2GYGsIZiWsDU0aNHsXjxYvTt2xdpaWkW8Ci4LijZqGUmEAcC4NzN7+HU2wNEAAi8NSb\/BxGJNWQDUgrAAKoOW4WKze+TLWfGDAxgM0ZNwzY7ATxkyBDcfLP89JCGVbMpAyqQnp6O4cOHg\/tDYMGZtvEovj5bzauRhxpF4aFGkYFVonPplOgiJJ09hKzZ\/ZA6bitiGrVRVWPORxOQs\/IlBrAP1RjAqrqU9TIfO3ZM\/MH97rvvrOcce8QKhFiBiw0642KD0tOn4edPoeL2V0H\/a+TUJzUHfaqfQflGbVB14EJEJqubEaFLNbJmXr7PmUfAniPNADbyJyBIbSMI0z9OrAArwAo4FaARcLVyRYiqmqYavmRDugnLbwDzGjB3SFaAFWAFWAFWQJ0CvAtaXi8eActrxDlYAVaAFWAFVCpQajNVTDxqvPA5ouu1kLUiBXdYfDJqvLRZvKBkxcQAtmJU2SdWgBVgBUKsQNGpP3BsdEuU5Fxe3lK6m1l6E1ZUagOkTvgakZWrhtgbfapnAOujK1tlBVgBVsDWCrjfBZ3Q9UUk9Jggq8nJN3sj\/5sPRL4KN3RCyvDVln2SkAEs2x04AyvACrACrIA\/CkgfVih\/7e1IHbkOYdEVvJri15D8UZnLsAKsACvACrACbgpI3wMOi45FtefXosI\/bveqE78HzF2IFWAFWAFWgBXQQIGS87nImHwXLh38VliLSKmLGmM2IjKlThnrRZmHkT65PYozD12efr6pK1KeWW7Z6WfykaegNehkbIIVYAVYAVbAswLnv9+AjFe6Iqz4ksgQlVIPVYd8iOirmgJhYYDDgQt7v8CJuY+iJPtPkScsNgmpL3yG8nWbWVpWBrClw8vOsQKsACsQYgUcDmSvnoQzK8cBjhJXY8JjKgMxccCFPJRc+PsRC0dENKoNXIiKrS\/fomXlxAC2cnTZN1aAFWAFjKBASTFy1s3AmRXj4Cg877VF4RUTkfzEO4i9qcvl0bHFEwPY4gFm91gBVoAVMIoCxdnpyPn36zi3cw2KT\/8OlBSDRrxRqVchvvVDiO8wULwBbJfEALZLpNlPVoAVYAVYAUMpwAA2VDi4MawAK8AKsAJ2UYABbJdIs5+sACvACrAChlKAAWyocHBjWAFWgBVgBeyiAAPYLpFmP1kBVoAVYAUMpQAD2FDh4MawAqwAK8AK2EUBBrBdIv2XnyUlJdi3bx8+\/PBD7NixA+np6eI3tWvXRps2bdCtWzfUq1cPYTY4g2eV0Ofk5GDIkCHYvn07pkyZggcffFCVa3r0CT1sqnLKJpmzs7OxefNmbNiwAT\/99BNOnz4tPK9RowYaNWokPs8tWrRATEyMYkX0iJ0eNhU7ZOCMDGADB0frptEf6kmTJmHNmjVeTUdGRqJXr14YOnQo4uPjtW4C29NYAYfDgffffx\/jx48XltUCWI8+oYdNjWUzvblLly5h8eLFmD59OgoLC336c+WVV2LcuHG49dZbZb9Y6xE7PWyaPoB\/OcAAtkokZfw4efIkBg8ejF27drlyxsXFIS0tDQUFBfjtt99QVFTk+l2XLl3EH\/WKFSvaRCHzuUnw3bZtm4hrXl6eagDr0Sf0sGm+yOjb4vz8fPFFeuXKlaUqSkpKQs2aNcXn+dChQ6XATF+sJ0+ejPvvv98rhPWInR429VU3uNYZwMHVOyS1EVhffvllLFq0SNRP4B09ejTuvfdeREdHi5\/RH3D6\/axZs1wgpmnNgQMHgj68nIylAMH3yy+\/xPDhw13wVTMC1qNP6GHTWKqHvjXuMx7UorvvvhvPPfccatWq5YIrxYK+nE2cOBG\/\/\/6763M\/b948NG\/evIwjesROD5uhj4C2LWAAa6unIa3t2bMHDz\/8MM6dO4fY2Fi89dZbuO2228q0lT7cq1atwpgxYwSEk5OTBZQbNmxoSL\/s2iiKzYoVKzBhwoQy049Kp6D16BN62LRrjL35ffToUTz22GNihEuJviAPGjTI9UXavdyRI0fw1FNP4cCBA+JXnTp1wtSpU8usCesROz1sWq0\/MICtFlE3f2g6auzYsa7pqh49eoj\/do583d2\/cOECRo4ciXXr1olfPfPMM+IDzpuyjNFRsrKyxB9Qb+v4SgCsR5\/Qw6YxFDdWK+hzSTNTlOiL8fz588WGK19JWqZ69epYuHAhrr76alcRPWKnh01jRUKb1jCAtdHRsFb+\/PNP9OvXD4cPHxZTye+8847H0a\/UgS1btuCRRx4RP7r++uvFh7xKlSqG9dEODaPZC4LuzJkzXTtdo6KiMGzYMHz\/\/ff44osvhAxKAKxHn9DDph3iqsZHmvmgKeUlS5aIYr1798aLL74ou0QkjQ2VW7p0aalpaD1ip4dNNVqZJS8D2CyR8rOddDSlT58+onSDBg2wYMECpKam+rQm\/fDQlDXtsm3SpImfLeBiWihAU86jRo1ymbruuuvEzta6deuK9XznjIUSAOvRJ\/SwqYVuVrJBo0r6MkzHjk6dOiW+WPfv31\/WRTqqNGDAAOzdu1fkpZ3TtMnSmfSInR42ZR01YQYGsAmDpqbJ7733ntgxSaljx45i+pKg6ivl5uaKaWf6EHn6wKqpn\/Nqo4ATwLTT9dFHH8VDDz0k4ui+ZKAEwHr0CT1saqMcW6F1YNoDQl+sKbmPgPWInR42rRhJBrAVo\/qXT7Spiv4gv\/vuu+InSqes3P+o00YPmurkFDoF1q9fLy5NoWMkiYmJroaoBbAefUIPm6FT2no1+1oD1iN2eti0XlQue8QAtmpkgTKjo8cffxwjRoxQ5PG0adNARxYoqSmnyDhn0kwBtQB2z68mtt76hB42NRPI5oZo0x5t2tq5c6dQgr7A0TqycxOmHrHTw6ZVw8gAtmpkPQCYzgJ3795dkcfSP7a0XkRlve2cVmSQM+miQKAA1qJPuLdBC5u6iGUzo7Rpi8710z9KtGQxd+5ctGrVyusMihax4\/6gvOmRAT0AABAqSURBVKMxgJVrZbqcdGPOs88+i02bNom2K1kfdDopBbC3s4OmE8SCDVYLYD36hB42LRiqoLrkfqafKn\/66afFP+nFOnrETg+bQRUviJUxgIModrCrct\/9yAAOdgT0r08tgPXoE3rY1F8569ZA8KUja7Rr3nlPdOvWrfHaa6+Jy3WkSY\/Y6WHTqtFiAFs1soC4+You1aANPDwCtmag1QJYjz6hh01rRkt\/r2jamXY508kH593uN910E958801UrVq1TAP0iJ0eNvVXLjQ1MIBDo3tQalX7x1naKF4DDkqIAq5EbYzV5lfSJ\/SwGbAwNjRALyTRGu+cOXNc8KX1Xpr58nZblh6x08OmVcPJALZqZHkXtIUj+7drav\/Y6bFDVQ+btgiehk56evKvbdu24gUk92lnabV6xE4PmxpKZShTDGBDhUP7xkhHsnS9JK0Lyd3r7P4Boo1cdDEHJ+MpoBbA5IEefUIPm8ZT25gtomtm6Zy+86YramXfvn3Fz5Q8J6pH7PSwaUz1A2sVAzgw\/Qxfevny5eKqQkpKdzO734RFxxjoFi1OxlPAHwDr0Sf0sGk8tY3Xou+++06A9vjx46JxtMOZvjDTNZVKjw3qETs9bBpP\/cBbxAAOXENDW5Deydq4cWNxK5avKSlyRnp1HV19SK+nXHvttYb2066N8wfAevQJPWzaNaZK\/KadzvTe7+DBg13vQdM736+++iratWsnO8slrUOP2OlhU4kuZsvDADZbxFS21\/1hBXrft2nTpj6t8GtIKkUOYXZ\/AKxHn9DDZghlNXzVBDi6IjYvL0+0tXbt2uKYkT+PpugROz1sGj4ofjSQAeyHaGYq4v4uJ13KTlPS0sP4Un\/4PWAzRbfsdaNKznrr0Sf0sGmuSASvtQcPHsSTTz4pZqoo0RfqGTNmoFatWn41Qo\/Y6WHTL+cMXogBbPAAadE8+rZMH1g6n0fgpZ2RdCes+2Ys99tzaKqaRsz08DcnYyrgzwiYPNGjT+hh05iqh65V7nc708iXjh7Vq1cvoEbpETs9bAbkpAELM4ANGBStm0TnAydMmIBly5YJ0\/SQOz3kfd999yEmJkb8jOBMT4jRhivnAX66xJ2mubyNlrVuJ9tTr4C\/ANajT+hhU70i1i1BX5Dpbe7x48cLJ2mT1SuvvIKWLVuqcpo+887PvbOgHrHTw6YqR02QmQFsgiBp0cSTJ0+KDRu7du1ymSMQ04PulA4dOuS6to7+mx5goA+6kmMMWrSPbfingL8Aptr06BN62PRPGeuVyszMFC+T7du3LyDnvC1T6BE7PWwG5LzBCjOADRYQPZtDH2C6om7Dhg1eq6HRbq9evTB06FDEx8fr2Ry2rYECgQCYqtejT+hhUwOpTG9CurM4EGd87RPQI3Z62AzEfyOVZQAbKRpBaEtJSQl2796N1atXY8eOHeKRd0q0ltSmTRt069ZNrCfJXdYRhKZyFQoUCBTAVIUefUIPmwrksHSWFStWiIt0Ak1yG\/X0iJ0eNgPVwQjlGcBGiAK3gRVgBVgBVsB2CjCAbRdydpgVYAVYAVbACAowgI0QBW4DK8AKsAKsgO0UYADbLuTsMCvACrACrIARFGAAGyEK3AZWgBVgBVgB2ynAALZdyNlhVoAVYAVYASMowAA2QhS4DawAK8AKsAK2U4ABbLuQs8OsACvACrACRlCAAWyEKHAbWAFWgBVgBWynAAPYdiFnh1kBVoAVYAWMoAAD2AhR4DawAqwAK8AK2E4BBrDtQs4OswKsACvAChhBAQawEaLAbWAFWAFWgBWwnQIMYNuFnB1mBVgBVoAVMIICDGAjRIHbwAqwAqwAK2A7BRjAtgs5O8wKsAKsACtgBAUYwEaIAreBFWAFWAFWwHYKMIBtF3J2mBVgBVgBVsAICjCAjRAFbgMrwAqwAqyA7RRgANsu5OZweNq0aZg3b54uja1evToWLlyImjVrYuTIkVi3bp2op1OnTpg6dSpiYmJ0qddKRi9cuFBKO0++LV26FM2bNzec29K+dd1112HBggVITEzUvZ1m1kx3cWxaAQPYpoE3utsMYGNHyMwwYQAbu2\/ZqXUMYDtF20S+MoCNHSwGsPr4mFkz9d5yCSUKMICVqMR5gq7AgQMHkJWV5bPegoICvPPOO\/juu+9EvqSkJDz33HOoUaOGz3KRkZG45pprEBUVxVPQfkbWHSZmmr4P1QjYk9QrVqzAqFGjXL8y6rS9n92Ei8kowADmLmJaBdwh4Fzbvfrqq03rk1kabmYAG0ljBrCRohH8tjCAg68516iRAgxgjYT0wwwD2A\/RPBRhAGujo1mtMIDNGjluNxjAoesEDGBttGcAa6OjWa0wgM0aOW43AziEfYABrI34DGBtdDSrFQawWSPH7Q4YwEohIt20M2XKFDz44INC\/WPHjmHVqlXYvHkzaNNYUVGR2AhGZ1979eqFZs2aITw8vFSknGW++OIL\/PLLL+J3tGmsdevW6Nu3L+rVq4ewsDDF0c3Ly8O3336LTz75BD\/\/\/DPS09NFWWpH06ZNce+99wrbsbGxim0qyahUO6etnTt3omfPnuI\/pRu2zp07h40bN+LTTz\/Fjz\/+CPKHNsddf\/31uP\/++3HXXXeVabuzDPlMG\/AKCwsRFxcnyvTu3Vv4Gx0d7dUNJZuwfv31V\/Tr1w\/Hjx+HdG9BSUkJdu\/ejdWrV2PHjh0uvSmGLVq0QNeuXXHjjTeCNvopSQxgJSpZNw8D2LqxtbxngU5BK4WIO4AJaosXL8b06dPFH39vqUuXLnjhhReQkJCAS5cuyZahP9oEkGeeeQYVK1b0GT+yt3btWrzyyis4ffq0z7wE4xEjRggYKwWDXOdRqp03AL\/88sviCwO16\/fff\/daXd26dfHqq6+CLsxwOBwCfnJl6IsPfVGqU6eOR7v+Arh8+fKibvoy4SvJ1S8tywCW62nW\/j0D2NrxtbR3oQDwSy+9JEY97777riJte\/TogWHDhmHWrFkCwErSkCFDMHDgQK+wzM3NxcSJE8UoTE2iETa1RQ7uSmwGAuCOHTvijjvuAGlJI165RBCeP38+9u3bJ77QKClDsxAzZ85EcnJyGfP+AJjA+8Ybb+DIkSNyzRW\/d7Y5LS3NZ34GsCI5LZuJAWzZ0FrfsVAAuHbt2vjzzz\/FdHPbtm0xePBgMW1crlw55OfnY\/369aDRnRMSNPVLx6JoepVGnw888AAeeeQR1KpVCxEREThx4gSWLVuGuXPnCpuUqA66HtHTH2+qg8C1Zs0aV4BphEvQ7tChg5h6pinsnJwcbN++HTNmzCg1wnz66adB\/wIdCQcC4KpVq6K4uFiM3K+88koxqrzlllsQHx8vNNi7dy8mTZok\/teZaAR8+PBhoSuNMMkHmmKna0NJk2+++QaTJ08WU8bORACm6W73pBbApBXFimYdaHr8scceQ7du3cTSAS0xkNYUd6pPOhsxYMAA4ZsvrRnA1v875ctDBrC9429q70MBYKdg3kBG06QrV64sdbkClaE\/wgQIWtd0X+OldcU5c+YIWDoTTW\/TFLZ7cv+DTV8CJkyYgJSUFI+xJDiMHTtWAMLZjrfffhu33357QLEPBMDOim+66Sa8+eabICC7J1pTJ4BJgUp5SJPx48d7HMXT9DR9uXF++aG89GXIfT1YLYCdbaMvRm+99RYaNmzoUbuDBw\/iySefdI2SGzRoIL5IpaametWaARxQNzR9YQaw6UNoXwdCBWCCHgGSNv54SjRCpg08NGJzpocffhijR4\/2OhqSbvqhMo8\/\/rgYPUlTZmam+DlNxVIiENDUrNzNX1SOvjDs2bNHlPvXv\/6F1157LaCp6EABTFPDBKdrr73Wo4Z0yxl9caAvM84k5697m7w9tOAPgOkLlNwXF\/ryRYB+\/fXXXV92PvzwQzFi95YYwPb9+0WeM4DtHX9Tex8qAHub2nSKSbt06ZUl56iTpqEXLVokpky9pezsbDHic067ehq9kT0CqTPJtUNal\/QPPbXn\/fffR5MmTfyOf6AAppkAWsf2tVt5+fLl4kuLM9HmtEGDBvncJS6FK03zk+40cpUmfwDcqlUrAVeaJveVaNq\/T58+riy0V8DXbAMD2O8uaImCDGBLhNGeToQCwDSdSE8Z0rqvt+TeLiVTke7Qdr9bmdZGCVhLliwR1aq9dpOmR2lUnpGRIcrTZqb+\/fv73XECBTBNDXfv3t1n\/e5wkoMZGZNC25tG\/gCYprbpzma5I2I\/\/PCDADDFk5L02JonZxnAfndBSxRkAFsijPZ0IhQAbty4sdgB7Wl3rTMK7u2ikS9NFdNxJKXQdgcw7Xym0R+NsCgpsSmty32ETZuIaO2YNo\/5kwIFMI3AaVTpK0nhRFPActO5ZEtaRksAy4HU6Yf7UoJcOQawP73POmUYwNaJpe08CQWAlTzgrhZOFDi5MvQyFI3CnOu\/gQY70NeL5Nrr3j7pRRz0OyWv\/iiBqXs9Ssr4MwKWAykDONAeac\/yDGB7xt0SXtsJwO4jq0ADqOSLhK86GMCe1eERcKA9017lGcD2irelvGUA+x9OBvA8IZ43HdSClEfA\/vdFO5dkANs5+ib33U4Adj\/aRFdWvvjiiwFfqOFvF+ARMI+A\/e07XO5vBRjA3BtMq4CdAEwXatANTM6zvHQWmS7u0OJaSX86AAOYAexPv+EypRVgAHOPMK0CdgKw+8UUvq6rDEZAGcAM4GD0M6vXwQC2eoQt7J+dAExhdD+yMm7cuFKXPvgK9f79+\/HUU0+JiyToLHP79u1x9913i7uN\/UkMYAawP\/2Gy\/AImPuARRSwG4CPHj0qbstyvshDo2B6xMHXpSAUanpEgM780qMPzqTmFi1P3YUBzAC2yJ+RkLrBI+CQys+VB6KA3QBMdw3TfcT0Pq4z0YUctBZM1y56SnSDFp25pdeFnK8t+XqqT2k8GMAMYKV9hfN5V4ABzL3DtArYDcAUKNqMNXz4cGzZssUVN3qCkEbGnTt3Fq8i0RN5pA3dK\/3OO+9g69atrrz0gMTs2bNlb6GS6xQMYAawXB\/h38srwACW14hzGFQBOwKYQkG3Yo0ZMwabNm1SFRmCL42Eae1X7k5jOcMMYAawXB\/h38srwACW14hzGFQBuwKYwkG+00s\/NCXtfP\/WV5joSTx6gIHustYiMYAZwFr0I7vbYADbvQeY2H87A9gZNoLvt99+i08++QQ\/\/\/wz0tPTxa9od\/M111wDWu+lqen69euLqWmtEgOYAaxVX7KzHQawnaPPvrMCfiqgFsB+VmP5YvwakuVD7NNBBrC948\/eswJ+KcAA9ku2MoUYwNroaFYrDGCzRo7bzQqEUAEGsDbiM4C10dGsVhjAZo0ct5sVCKECDGBtxGcAa6OjWa0wgM0aOW43KxBCBRjA2ojPANZGR7NaYQCbNXLcblYghAq4A9hTU+gGLtqFzemyAqwZ9wR3BRjA3CdYAVZAtQIME9WSMYDVS2b5Egxgy4eYHWQFtFeAAaxeU9ZMvWZWL8EAtnqE2T9WgBVgBVgBQyrw\/0m+4kA0SfQ6AAAAAElFTkSuQmCC","height":193,"width":320}}
%---
%[output:78729cd3]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAeAAAAEiCAYAAAAh9AEqAAAAAXNSR0IArs4c6QAAIABJREFUeF7tfQvYTlX6\/sPIoVCfnA\/llBhKKpP5RyExMw1FdCA5VBTp4BCiHFLIoatCSCilosNEaaYUDfMrqqFJqVQUIceiRmHyv+5V+7W\/7X2\/d+1373fvvfZ7r+v6Lodv7XW477XXvddaz\/OsQkeOHDkiTESACBABIkAEiECgCBSiAAeKNysjAkSACBABIqAQoABzIBABIkAEiAARCAEBCnAIoLNKIkAEiAARIAIUYI4BIkAEiAARIAIhIEABDgF0VkkEiAARIAJEgALMMUAEiAARIAJEIAQEKMAhgM4qiQARIAJEgAhQgDkGiAARIAJEgAiEgAAFOATQWSURIAJEgAgQAQowxwARIAJEgAgQgRAQoACHADqrJAJEgAgQASJAAeYYIAJEgAgQASIQAgIU4BBAZ5VEgAgQASJABCjAHANEgAgQASJABEJAgAIcAuiskggQASJABIgABZhjgAgQASJABIhACAhQgEMAnVUSASJABIgAEaAAcwwQASJABIgAEQgBAQpwCKCzSiJABIgAESACFGCNMbBlyxbBDxMRIAJEgAhEC4GqVasKfkxMFOA0rEF4z5y8QeU66W\/Xmcgx20wEiAARiC0C5513nkyYMMFIEaYApxmW77zzjvxlwQGVa8bFhaRaXvHYDuRc7NiqVavkwQcfVC9wlSpVchGCWPeZ\/MaaXrH4nT9\/vjRp0sS4zlKAXQjw2uF\/lFPKhCvAh3dukiLlqsv\/9u2Ug5s\/khL1mxs36KLUYHxgde7cWUx9gaOEZRTbQn6jyIp\/bTKdXwqwYQL8ZadC+VoMMa48apkSZSb3CGzatEkef\/xx6datm1SvTgzdIxjtJ8hvtPnx2joKsFcEI\/68fQs6zBXw\/uVzZefUHinRqrnwSMSRjGbzfvrpJ9m2bZtUqlRJihcPd3cjmgiZ3SryazZ\/6VpPAU6HkOG\/twvwoj6NpGntkwLvEbadt45oIfizoFRp5DJuSbtkhxO0S8AMy05+DSPMZXMpwC4BMy17FAR4z5N3yHcvTUhAh+1mnP1iVWxP+P9Tpm00DeJQ28sJOjj4w3DnO3z4sOzatUtOOukk7nAEQHXQLkEU4ABIDbMKuwBPvbqeXN24YqDNwar36z418tVZru8cKdW8uxz4aLlsG9ki3+\/yOo2UvCtGBNpGkyujAAfDHsR30KBBymqVKb4IBO0SRAGO71hSPQtbgPcuGCV7F47Mt\/q1r3L3PnOX7H1+TD4W8HsaZekNTAqwHk5ec1kTJd29vCIZ3efDcAmiAEd3PMgPP\/wgI0aMkBdffFF69+4tgwcPdt3aMAW4oNWvvSNOy2isjrFKZkqPAAU4PUZ+5DB9ovQDg7iXEQbHYdTpJ4+xdUM6cuSILFiwQIYOHarwMlGAnVvMqc54k1lI0ypa7zWhAOvh5DWX6ROl1\/7nwvNhcBxGnX5yGVsBXrdunVx33XWyc+dO3wR4cJvqMrhN\/vNYP8lwlrV1ZAv56aPlif8uaGWLc2K7lbR1TpzN9sWhbApwMCyaPlEGg5LZtYTBcRh1+slSLAV47969cuutt8rKlSsTWPmxAg5SgJNtPxfkZuQ8K0bHuQpO\/6pQgNNj5EcO0ydKPzCIexlhcBxGnX7yGDsBhtvB1KlTVXzfhg0bqhXw1q1bfdmCDlKAdbef7YPBeRZMv+D0rwoFOD1GfuQwfaL0A4O4lxEGx2HU6SePsRNgrHr79u0rZcuWlXvuuUcF2f\/ggw98EWC4IMEVKYjkXNHqGFYhUpbdN1jnmSD6EuU6KMDBsGP6RBkMSmbXEgbHYdTpJ0uxEuDt27dLv379lODee++90qJFC7nhhhuMFGDn+a\/OajaZXzBdkgp+XSjAfk4nqcsyfaIMBiWzawmD4zDq9JOl2Agwtp7vu+8+mTt3rnTo0EGtfg8cOKAMsUxbATvPf91EuKIxlrvXgwLsDq9Mc5s+UWba71x6LgyOw6jTT05jI8BLliyR2267Td1oM3PmTPXnnj17fBXg82udJIv7NvIT\/6RlOVeyxes3l8ojl2nV63RJcvOsVgUxy0QBDoZQ0yfKYFAyu5YwOA6jTj9ZioUAf\/bZZ3LTTTfJ5s2bZfr06dKyZUuFkakC7Dz\/dRte0mmMxW3o1K8MBdjP6YRb0MGgGc1awhDDMOr0E33jBdge7ap79+5y5513SpEiRbIiwNWL\/yAv9DpDKlbMbjzoPeP+nM\/\/N2\/Iq65uOdox5mI5\/OlRF6wSTbtIXu9Zfo6b2JQFAYbtQIUKFaREiRKx6VfUOvLee+\/JtddeK\/Pnz5cmTZpErXlsjw8IWGI4e\/bsrHOMdxY\/33zzjYoxbuq4MlqA7dGuzjnnHHn44YfziaPfK+DC\/90l\/Sp9rCaSrKW9W+TI2AuOFp9XVQoN\/ae76t57Xo4sGJTvmUL3f+mujBzJ\/fPPPytXtXLlykmxYsVypNfBdxN2GP379zd2ogweMfNqtAR48uTJygU0m+mJJ56Qxx9\/PFEFBTibaKco24p2hVUMfH+bNm2aL6ffAnxKmeLyt2sqZXUFjPPfveP+nOgHDLDKT17vGt2tXfOv5tyuol1XaOgDMNT79ttvFafFixc3tBfRbzYm5549e1KAo09Vxi20BBjieO6552Zcjs6D1go4jAsgdNqnm8fYFbA92hVcj\/BjbT1bnc+GAK8d\/kddbDPK5\/X816qUPsF68PMMWA8nr7lMP6vz2v9ceD4MjsOo008ujRVgC\/hMwBg7dqxceeWVWo+inr8sOJDIu2dy\/vt3tQpxkcnp\/5tpTGde0KAHOgVYDyevuUyfKL32383zOFrDrgx2+LB1v3r1asFF92PGjJHChQvL4sWL5YEHHpCDBw8K5rJWrVqp4j\/88EMVeAhYwxVz+PDhUrJkSTdVe8obBsdh1OkJJMfDFOA0aAYtwE4\/Xp0AHKm6QGvo9K8KBTg9Rn7kMH2i9AMD3TIgwPv371dCOmzYMNm9e7cS2ksuuUQFGIKowvj02Weflfbt2ythXrFihcD46bzzzlN\/wrbhqaeeksaNG+tW6zlfGByHUadnoGwFGCvA+PrDICwoff\/993LLLbfIRx99JN26dVPb1EiwdtW1eA1SgL0E4EiGA7eh078qFOD0GPmRw\/SJ0g8M3JaxYcMG6dGjhxJTxDZYtGiRnH322fLXv\/5V5syZo4T3+uuvl2bNmslLL72kxBoGhXgG8yPyNGjQwG21GecPg+Mw6swYoCQPGivAOiD4fQaMOnEGDGOsbCQvATiStYfb0OlZogCnx8iPHAVNlGX66wWZ8aMdUSgD8eQRVz5devPNN5XAwm0LYXXz8vLk8ssvl0OHDsndd9+t7jvHrW9ffPGF2m4uX768vPLKK2qhAYPUKVOmSOnSpdNV49vvwxDDMOr0DTARoQCnQdO5As6mAPtlgGXvErehCyaYAuzndJK6LArwUWx0BRgCCpeeGjVqyIUXXiiDBw9WrnK7du2SXr16yfr166V27doyceJEOf300wVb1zgDRjCiG2+8UfnHFipUKFExtqn79Okjo0aNUmfEfqcwxDCMOv3EjQIcIQH2ywDL3iXnNnSmRl1+DroolUUBDoYNCrA7AcbxGvymly5dqkTWCq+LUtasWaNiEfz4448yYsQI9XcI7b59++Tmm29W96DPmjUrERHQqvnvf\/+73H777Wq1DCNUp9eI15EQhhiGUadXnOzPU4BdCvCiPo2kae2T\/OQgUZafBlhWoYwNzRVwVgary0IpwO4E2Dr\/xV3mdpFFKTC+Gjp0qJxxxhkyY8aMRFwC6xn4s+P8t1q1ai5Z8pY9DDEMo05vKOV\/mgIcEQH22wDL3i1uQ6cmmStgP6eT1GWZPlEGg9LRWqyzXGw\/P\/bYY+pyGSQYV1nnv1jNYsVrbTNbZ8Zt27aVcePGaRua+tW3MDgOo06\/8EI5sRZgP4ByngHrnt+4rdtpgOXmCsJ0dWVjaztdnab8ngIcDFOmT5TBoPRrLbhaFdepzps3T7p27Sp33XVXYrt427Zt6oY3XDyDiFONGh29nW3SpEkqIiC2mBF1bMuWLWp1jGfgCYIfbFMPGDAgK8ZZYXAcRp1+jgUKsMsVcLYEOJtbxc6ySzXvLjgLZhKhAAczCkyfKINB6ddaLCOrtWvXHnOWCxxx5gvLaLuVs3VmvHz5cuX\/i\/NdBOzAOTL8hnEujOf69u2rLKc\/+eQTdb789ttvS61atZTBFp579NFHpW7duvLII49ImTJlXHU7DI7DqNMVKGkyU4AjIsDZsIC2d43b0MmJpgD7OZ2kLsv0iTIYlH6txTKywg1dzrNcyzLaaeWMFTH8fxH7AJbOsJCGqOKSESScG2Mlbb+u9bXXXlPW0vAnxmr54osvlldffVX5FM+dO1fq1avnqtthcBxGna5AoQB7g8u5BT24TXUZ3KaGt0KTPO3cJvYSAStZ45wGXrSG\/hUlCrDvQzlpgaZPlMGg9GstiGQFUbziiitk9OjRUrRoUfX\/dsto3PyGyFhWQgCOgQMHyscff6y2n7HaPeGEE9SvrXNjrIIh6KeddlrCZemFF16Q888\/X6655hq18h0yZIh89913gu3ssmXLuup2GByHUacrUCjA3uAKSoCzYQFt7zmjYiUfBxRgb++H7tOmT5S6\/YxiPuvcGPGkIaylSpVKuCyBF4S37Nixo1oFYxX95z\/\/WW1Tu3VTCoPjMOr0k2NuQadB0ynAiGCDc2A\/UzYtoK12MioWBdjPMeu2LNMnSrf9jVJ+a0sbW9M33HCDElZc9ACxxVky4kzjnBhW1BBenP86r3bV6U8YHIdRpw4WunkowBEQ4GxaQNu7x23oY8nmClh3qvCWz\/SJ0lvvw30aW9q4PQnCCrcmbC3jTBjb3NaZsBVF66233lJBP0488US1hW2PpJWuF2FwHEad6XBw83sKcAQEOJsW0NyGLphgCrCb6SLzvKZPlJn3PNwnrfNfXGsIv2HEje7cubPyJf7888+V2FapUkVwwX3v3r2lTp06cs455yiRvuiiiyjAWaaPAuxSgM+vdZIs7nvU984PfrJtAc1t6NQsUYD9GMHpy6AAp8coGzngU4yV75NPPimdOnVSMaRhKY3t5wsuuEDuvPNOtSVtGXHt2LFDCXXLli15BpwNQhxlUoAjIMDZtoAuaBvab2vrAMasr1VQgH2FM2VhFOBgcA6zljA4DqNOPzGmAEdAgLNtAc1taK6A\/Zw0MinL9Ikykz7n2jNhcBxGnX7ySgF2KcC4CxhXEvqVnBbQKLfmwiN+FX9MObSGzg8JV8BZG2r5CjZ9ogwGJbNrCYPjMOr0kyUKcMQE2M8Y0Km6xqhYR5GhAPs5naQuy\/SJMhiUzK4lDI7DqNNPlijALgUY2fdMbuEbB04XpOL1m0vlkct8Kz9ZQbwjmAKc1QGWpHDTJ8qg8TKxvjA4DqNOP7mhAIcswEFZQNu7ycsZKMB+TiI6ZZk+Uer0MdfzhMFxGHX6yTMFOAMBxhkwzoL9SM7VaF6nkZJ3xQg\/ii6wDG5D\/woPt6CzPtRUBaZPlMGgZHYtYXAcRp1+skQB1hDgv87dLL8cfzQwuZ8CHKQLkr2rvCOYAuznRJKuLNMnynT94+\/D+cgyfVxRgEMWYKcL0inTNgoMsbKduA1NAc72GLOXb\/pEGSRWptYVBsdh1OknPxRgDQFuN22tHC57eiLnoj6NpGntkzzzEMQlDAU1ktvQ3IL2PIg1CzB9otTsZk5nC4PjMOr0k2QKcA4LMC9noAD7OZkUVJbpE2VQOJlcTxgch1GnnxxRgDUE+LKJS+XgKecncuI6QlxL6DWF4YJkbzO3oSnAXsew7vOmT5S6\/czlfGFwHEadfnJMAQ5RgMNwQXJ2N9e3oWkF7ed0kros0yfKYFAyu5YwOA6jTj9ZogBrCHCH0c\/IT3UvTeQc3Ka6DG5TwzMPYbkg2Rue69vQFGDPw1irANMnSq1O5nimMDgOo04\/aaYAhyjAYbkg2bvs\/Ago1by7lOs7x88xFumyKMDB0GP6RBkMSmbXEgbHYdTpJ0sUYA0B7jhkqvz37J6JnDj\/xTmw1xTkLUip2uo8h0a+oFyhvOLnx\/MUYD9QTF+G6RNl+h4yRxgch1Gnn0ynFWBc6PzRRx\/J\/v37\/aw3ZVmnnnqqVKtWLZC6dCoBwdkQ4KBvQSqor7m8DU0B1nkLvOcxfaL0jkD8SwiD4zDq9JPJtAJ84MABGTJkiCxevNjPelOWNX36dGndunUgdelUEpQAB3ELUqr+Oo3BcmkbmgKs8xZ4z2P6ROkdgfiXEAbHYdTpJ5MU4DRoguArbhklPzS9I5Hz\/FonyeK+jTzxELYLkrPxuWoNTQH2NIy1HzZ9otTuaA5nDIPjMOr0k2JtAT7zzDOlW7duUqRIET\/rT5S1Z88etdLu2LFj5FbA2RDgqK06c3UbmgKcldf5mEJNnyiDQenXWo4cOSLffvutrFu3Tj744ANZvXq1VK1aVcaMGSOFCxdWu5EPPPCAHDx4UMaOHSutWrVSz3344YcyYcIEdfFFhw4dZPjw4VKyZMnAmh4Gx2HU6SegWgI8YsQI6dq1q5xxxhl+1n1MWePHj5dGjRpFToCvuv4W2dd6fKK9uAkJFzJ4SVHwAba3P1eDclCAvYxi\/WdNnyj1e+o9JwQYNjfAbNiwYbJ7924ltJdcconce++9SlR\/+OEHefbZZ6V9+\/ZKmFesWCGzZ8+W8847T\/35888\/y1NPPSWNGzf23iDNEsLgOIw6NeHQykYBTgMTCHYKMB7ZM7mFFsCpMkXBBcnZtlzchqYAexrG2g+bPlFqd9THjBs2bJAePXooMZ05c6YsWrRIzj77bPnrX\/8qc+bMUcJ7\/fXXS7NmzeSll15SYr1z5071DFbHyNOgQQMfW1RwUWFwHEadfgJKAdYQ4M6dO8t3lz2WL2ccBTgXryikAPs5naQuq6CJErtBuZSK179QStRvnrbLb775phLYJk2aSIsWLSQvL08uv\/xyOXTokNx9992yYMECufXWW+WLL75Q283ly5eXV155Rfr16ydNmzaVKVOmSOnSpdPWgwwQ+aJFi0qhQoW08ifLFIYYhlFnxgAleTCtAMMN6ZNPPpFTTjlFm8xMG\/jpp5+qQYaBFJVkEey3AEdxtenchi5ev7lUHrksKlRkpR0U4KzAekyhBU2UznchmBaFVwsC3cDTIF2CgE6ePFlq1KghF154oQwePFiKFSsmu3btkl69esn69euldu3aMnHiRDn99NPV2THOgOFJcuONN8qgQYMKFNSVK1fKv\/71L1m+fLk69oNwWzY+X3\/9tTz00ENK\/GGXo5PCEMMw6tTBQjdPWgHWLSiu+SyCcQb8y\/FlE93EGTDOgjNJYV9DWFCbnZNhpZHLtL7WM8EhCs9QgINhgQJ8FGcdAcYZb\/\/+\/WXp0qVKZLEFXb36r\/eEr1mzRq699lr58ccfBfY5+DtWrvv27ZObb75ZIKyzZs2Sli1bpiUXeXv27KlEG\/n37t0rzzzzjDp\/xrkyzp6vvPLKtOUgQxhiGEadWmBoZkorwDhLeOONN9SWRqlSpTSLTZ8N2yt169aVypUrp88cYo5sCLDTBSlMH2AntLkWmpICHMzLRQF2J8DW+e\/WrVvziSxKgfHV0KFDlVHsjBkzpGLFX29ms54pXry4Ov\/VCWiEVTYEF\/lPO+20RCMtvijA2X0\/0gowAnFgWwNfSTCF9yNB1GHNh62NbFtWe22vNRDhB3y47OmJ4hb1aSRNa5+UUfFR8wG2d8K5DY3f1Vx4JKN+mvAQBTgYlngGfBTnUi26CT66C0rWWS62nx977LHE6hdzp3X+e\/vtt6sVr3Vua50Zt23bVsaNGyclSpQosA5rlY1MkyZNyrfAogAH815oCTD8c5s3by6XXnqp8kPzmj7++GMZMGCA3H\/\/\/cYIcJl2d8uXhU9NdN3LncBRc0Fy8plLPsEUYK9vs97zpm8V6vXSn1ywu7nnnntk3rx5yv3zrrvuSpzNbtu2Ta677jrZvHmzPPHEE8pt00oQ0alTpyqDLCyYtmzZolbH9tgNKPvLL79U58j4f2xzd+nSRW644YZ8+SjA\/nCZrhRtAYbzd5UqVaRs2aPnoOkKT\/Z7THiw2jvhhBPUADJlBeynAEfhGsKCuMulbWgKcCZvsftnKMD6mFlGVmvXrj3mLBc44swXxlF2K2drNQuDKvj\/QlwxZ0NgrWAcCNQBwypYUtevX1\/9Dga2jzzyiDpitCcKsD5fXnK6EmAvFTmfPfHEE3NWgKPoA2znJ5duSKIA+\/lWpy6LAqyPs2VkVaFChWPOci3LaKeVM1bE8P\/9\/vvvpU+fPspCGlbQ5cqVUxXD2Ao7mdi+vvjiiwXjHv\/+\/PPPlYEXFlcUYH2O\/MpJAU6DpDVxNOx2n7z1fYVE7sFtqsvgNjUy4iEK1xCma7jzIyGv00jJu2JEuseM+z0FOBjKKMD6OCOSFYJsXHHFFTJ69Gjln4tkt4x++OGHVWQsKyEAx8CBAwXHe9h+xioZu4xIn332mdx0003Kyhl57K5MderUkVGjRqn\/owDrc+RXTtcCjNBmiMZiDQq3DcGXGqyqcZZs0hZ029sny7yNJya66+VOYKerTxSNnHLFGIsC7PYNziw\/BTgz3Lw+hQAbEPG33npLGXPBXxjp\/fffl+7du6vz4mRuRtyC9oq83vNpBRjFIC4pwqBNmzZNYAQAErG9gXODTIQYxgEwwgL5ppwB+yXAUfYBdg6ZXPAJpgDrTRRec1GAvSKY2fPW1nTNmjUTls4wxLrvvvtU1CyIcrJwlRTgzPB2+5SWAFuF4msKztmPPvqovPvuu3LyySerUGiwjnbjI4wBAAtomMvnmgBH2QXJOXhywRiLAux2ysgsPwU4M9y8PmXhbremxnlw37595ayzzhJcgAPXJvgM28NQUoC9Iq\/3vCsBtor85ZdfBGEjcQ6BSC0gDqbsMI9HYA2deKJ4DitpHWdxva5kJ5c1EPveM03uff\/oOUmmdwKbFO4xF4yxKMDZeW+cpVKAg8HZWQuuNIRxVsOGDdUKGEZXuMoQ1xy2a9dO6tWrp1yVLrroonzzNkQa58gMxJFd3jISYHuTtm\/frszecZ6LrWrcTQkHcYirHz7D2e1++tKtiWP81NnS+\/WjASkyvZIw6j7ATkScxliIYYtQenFJFOBgmKQAB4OzsxbsNsLNCNbTiJAFf1\/EdEAkLfzujjvuUHGmrbkagZdwPozIWMuWLVOXQOAZCHi6wB5hcBxGnX4y6VmArcYgLumSJUvUOfFXX32l\/MxALu6nzOSc2M9OeikrlQCjzExuRHJu6+rEhfXSfq\/Pxt0YiwLsdYToPW\/6RKnXy9zOFQbHYdTpJ8u+CbDVKHxVrV69Wn1xAZxKlSopvzRsd7g5J\/azk17KshP8lwUH8hWViQBH3Qc4GVZxjoxFAfbydug\/a\/pEqd\/T3M0ZBsdh1Oknw74LsNU4XI21adMmFSz8xRdfVP8N\/7Ru3boJHMx1zon97GimZfktwE4xO2XaxrRxYTNtu1\/PObfNUW4UXacy6S8FOBPU3D9j+kTpvse590QYHIdRp5\/MZk2ArUbC5Qjb0rhxAwkh0nBODIdwXK8V9XNiO8E3Lj0iX+\/5KYG\/2ysJnS5IJgmZ0yUp6lvnui8JBVgXKW\/5TJ8ovfU+N54Og+Mw6vSTzawIMFa\/iL4CdyX4D2Nb2pngQ4wbkUyxgp4\/f744BdjtjUgm+QA7+XKeXRev31wqj1zm51gMpSwKcDCwmz5RBoOS2bWEwXEYdfrJkq8CDPck+Acjtigs6JKlyy67TFnVmWIlbScYbkj\/+uK7RLfc3ohkkg+wk7tkq\/dKI5dJifrN\/RyPgZdFAQ4GctMnymBQMruWMDgOo04\/WfJFgK0AHbgKC\/5lzoRtZ5z9whnc6fDtZ2eyUZadYISifPrd7RkLsEk+wMmw3D6urfz3\/ZcTv4rDKpgCnI235tgyTZ8og0HJ7FrC4DiMOv1kyZMAw\/Xob3\/7WyJEpbNhsHrGrR0dO3ZM3MrhZ+ODKMtPATbNB9iJb7LAHKavginAQbxFojwiOnfuLDjKwVV6TPFDIAyOw6jTT+YyEuA9e\/aoFwlnvAi+4Uynnnqqcj3605\/+ZKTrkb0\/doJxG9L4f2xK\/NrtjUhOFyQTDZm23\/dn+e+av8dmFUwB9nM6SV2W9R4hdC1iAzDFD4FvvvlG3REQ5EdWTgkwLJohurBoPnTo0DEjCNFScO1Vs2bN0kZNwdVaCNAR9SAddoI3\/q669H16faLfbm9EMtEHOO6rYApwMEKAuQOT86pVq4KpkLWEggA+riZMmCBVq1YNpP7YC7Bl0fzggw+quM\/JLJqtcGXnnnuucjNKl1Dm5MmTpXXr1sZcxoCvOq8CbKIPcDIut49vJ\/99b3EsVsEU4HRvq3+\/hwjjJ8gEfr\/77jspW7as1twUZNviWBeENyjxBX6xF2DEBh0yZIgsXnx0wkXHIbSdOnVSAbtPO+00V\/68e\/fuVV\/D2I4y5TYkCPDhsnWl3bQ1iffGzYUMJvsA66yCTdxOR78owHGUgaN9Ir\/x5jcnBRgXOF911VVSpUoVV8L7v\/\/9T7744guZPXu2clfCBQ4mCXDlOmfJWWPeToxoNxcymOwDnOwVdvoFI4+J0bE4Qcd7gia\/8eY3JwXYD0pPPPFE4wQY1ptl+uf3b9aNB22yD3AyvpOt6PM6jZS8K0b4MTwCK4MTdGBQh1IR+Q0F9sAqzRkBxrYxtoxr1KiRMbimr4C9CLDTBzgO1\/olixFtmlsSJ+iMX2cjHiS\/RtCUcSNzQoBHjBghHTp08M1\/LxtnwFYUrieffFIdzO\/evVudUzdo0EDatm0rl156qZQpU8Y10U6CsQWdSTxo032AUwHnNCwzLTgHJ2jXr4RRD5Bfo+hy3dicEGCYleMmI7+s2\/y2gt65c6eMGjXZxnSNAAAgAElEQVRK3UecKlWuXFkmTZrk2gfRLwGOgw9wMmyT3Rds0iqYE7TrOc+oB8ivUXS5bmzsBRhuR7hYAeJbunRp1wCleuDTTz+Vk08+WbkHeEnwJ8YK3bryELGmu3TporbK9+3bJy+\/\/LI8\/vjjakWMyFxz586VRo0aaVfpJLjt1DUZxYOOgw9wKtCcfUM+E65ZRDs5QWu\/CkZmJL9G0qbd6NgLsDYSLjOuWbNGia\/X25Cw6r3tttuUf3K\/fv3Uj9MXGXXh\/7du3ap8jydOnCglS5bUarGTYATiyCQedFx8gJOBlswg6\/iz\/yIVh76ihXGYmThBh4l+9usmv9nHOMwaKMAZoo+IWjBq8uKGZPdRRlm4DCIvL++YFmHLe8qUKfLAAw+oVfecOXPU2bBOchI8\/h8bXYejjJsLUjLcTDXI4gSt8xaYm4f8msudTsspwL+hBCMoRJz56quvBJc0FJQ2bdqkXJBwJutFgDdv3iw9evSQL7\/8Uvr37y8333xzymotopDBTaxSJ8FY\/boNR5kLAgxck21FR\/08mBO0zjRnbh7yay53Oi3PeQHGKhSCNn36dHXOqpv88APesGGDOv+FEEN8ESAkVQpTgOPmA5wK42Rb0UXKVVfnwVFNnKCjyow\/7SK\/\/uAY1VJyWoBx7optX8SJdpv8EGA3dSL61pgxY+SEE05Qq29dQywnwSs\/\/851OErn9mwcfIBTYZ\/sysJSrW6Qcr1nuqErsLycoAODOpSKyG8osAdWaU4L8LZt2+S6666TTz75RCpVqiTXXHON1pby2rVr1Y1KjzzyiFZ+r2xu375dGWG9\/\/776tw51VlxsnqcBMMH2G04yrj6AKfi5fslD8ruObfl+3VUt6I5QXt9u6L9PPmNNj9eW5fTAowtYJzBwvBpxowZAl9bnQSjKKyaL7rooqwLMFbpDz\/8sPpBGjt2bIFb1c72JyPYbTjKOLsgpeLblPNgTtA6b6y5ecivudzptDynBXjXrl3Sq1cvZVF89913u7ruC5c3wy\/XT99iJ2EQ+ueee06GDRum3JTat2+vAnbouiChPOdF4hUrVpSzH8p\/prn1vj8WOFb2jPuz\/PTR8kSevCGvSon6zXXGl7F5Du3YJDsH1Dum\/eWGvy7Hnd40Mv3CBI0dkgoVKqS9wzoyjWZDtBEgv9pQGZUR7yx+oCO4Wc+NYW2UOlroCFQqw2S598Dy+d5775VixYpplwTDqRIlSngOxJGqQrQNwTmGDh0qhw4dkj\/84Q\/y0EMPSfny5bXbaBdg66Fu3brJwiMtZeu+w4lyFnerKpVLp74H+cgdNfPVWWjoP0Xygrmw2lVn\/c78xTtyZEbnY0ot1Hu+SK0mfteWUXk\/\/\/yzIJJauXLlXI3fjCrjQ4EjQH4DhzyQCmHHgwBLVspJAUbnd+zYIcOHD5fevXvLOeecow3+U089JWeeeWZWtqDhErVw4UK566671Mq3YcOGyuWpZs38QqjTWGsFjHCcuH4RK+CbXtopq776IfH4C73OkCbVkwf2cK4EYRVcfvJ6napjkSfVSjgquwCw4v\/2228Vr8WLF48F5uzEUQTIbzxHg7UCXrVqlTrOzFkBBr0IVYkt6JtuuimtoOJGpA8++EAFxRg3blza\/G6HD754Ydw1bdo0Jb74KJg8eXLGEbeSnTG4iYbltAqOuluOW7x183\/ZqdAxWaMQrpJnhLoMmpmP\/JrJm26rc\/oMGCBhq3f16tVK5N59911d3CQbbkiI\/XzPPffI888\/r9rRqlUrtTWO7cVMUzKC3UTDyhUf4HT47l82V3ZO6xE5EeYEnY45s39Pfs3mL13rc1qAnUZO6cCy\/95vAcZW+MCBA2XlypWqmssvv1xtQXs18kpGsJtoWLnmglTQGEglwuX6zhH4RoeROEGHgXpwdZLf4LAOo6acFmDLChp+vUh169ZVfrawbi4oIXQkgENwDC+hKK06IL633HKLWokj9e3bV0XGcmMUlqq9XgV459Qegiv7rJTXaaTkXTEijLEaiTqTBepAw8K6R5gTdCSGRdYaQX6zBm0kCs5pAbb8gGFljFCUZ511lhQqdOxZn5MpnM3ef\/\/90rZtW88CbL+OELcgwSQdvsnOG5EyHS3JCHZGwzqlTHFZOzy5K1Iu+gCnwzqVCON8vPKoZYI\/g0qcoINCOpx6yG84uAdVa04L8MaNG6V79+7StGlTFZO5aNGi2ri\/9tprUr16dalTp472M86M2AKH8MNCGSnVdYQZV2DzA3Za2ekG43BeQxjViFBeMMrk2WRxo8PYJeAEnQl75jxDfs3hKpOW5rQAw+gJW72409etAGcCtvMZawWOe34RChOrX9wxnC6deuqp2lbRqQh2CjBWwFgJO5PT+rfmwozdrtN1y7jfQ4R3TuspB9YtO6btWAXjbDjbAUs4QRs3bFw1mPy6gsu4zDktwGDr2WefVXGd4fYDEdRNfgTisC5Y0K3TyucmHGUqghEPGnGhrZRMgHPlGkK3+DvzJ7tL2MpTomFrKdd7Rta2pTlBe2Uv2s+T32jz47V1OS\/A1hnshRdeqM50dc6AAbrXQBwHDx5UvscLFixwzaEfAqzjC0wXJH1qUp0LWyUcf\/ZfpOz1U30XYk7Q+hyZmJP8msiafptzXoABFayQ4fLTuXNnLaMqrH4RhOPOO+\/Uyq9Ph\/85UxGsI8CwfoYVtJXCsvT1H5XslIgdg\/3LHpe9C0emrMDvrWlO0NnhMiqlkt+oMJGdduS0ANuvI3QLr99+wG7r182fimAdX2D6AOuinD8fhHjH1B75LrBwlgQhzut0t5Ro0MLTqpgTdGYcmfIU+TWFqczamdMCDHciRJ6aN2+ea\/RyQYCdLkhhBpxwTVAEHsC2NHYQIMgFJWtVfFz56q7FmBN0BIjOYhPIbxbBjUDROS3AwP\/NN9+U66+\/Xk4\/\/XRp3ry5livSunXr5D\/\/+Y9vgTiyOQ5SEazjC0wfYH+Y0VkRWzVBjI8\/t62ccF4H0RFkTtD+cBTVUshvVJnxp105L8A4z0VQjdGjR0teXp4Wqlg5wxDqsssuM\/YMGB1N5wvsdEGKwuUDWgRFNBOE+L\/\/XiK7ZvXVbiEEuUh5bFePSCrInKC1oTQyI\/k1kjbtRue8AOO6L2xBd+jQQcsH10LWj0Ac2ix5yFgQwQW5ItEFyQPoGo\/CwO2nj97KF+ZT4zGVxYq0dfy5f5XfndJQvtv7nZxUs4HknfMn3SKYzxAEKMCGEJVhM3NegDPEzZjH3Ajw1KvrydWNK6q+8RrCYCjGhw6w\/t++nbJn3h2+VGqtmotWrSfFap4rUrhwPuHG1rZdyH2plIVkBQEKcFZgjUyhsRdghHvcv3+\/uqzcTajJTBiCT\/Hvfvc7KVGiRCaPZ+WZggguyBWJPsBZoSNtoRBkCKjz\/D3tgx4y2GNXY7u7UOEiclzVulK46PFS6LhiUqRCzQJLTxX72hL6VA8HGTPbAzyhPkoBDhX+rFceewHGFjPCTHbt2jXr57Xjx4+XRo0aSevWrbNOnG4FBRFckCsSXZB0Ec5+PnwMHfpmvfzwfwsKdG3KfkuiUUOuCTdsTnA5Cz6OmOKHwKp33hHE2MdNfKalQkewxC0gUYDfUQFGnJcxALKCBJguSNF+FSDKWC3\/9N8fZf+K+VK8WFGKc7QpY+uIQEoEdgx4mwLsdXyYtgJGLGgYYtnTnskt1D\/pguR1NATzvHOL0vI5PrRjkxLowzu+Ug3B3w\/95o98+LffBdNC1kIEiEA6BGItwMOHD5dLLrlEWrRooR3rOR1gzt\/jnHnAgAHSsWNHY7ag0YdUtyLRBcntCAgnvx9nhPZAIRBuS7CtHlkinvi3I7CIJexOBCD0BaV0AUrCQZS1EoHgEYi1AA8ZMkReffVVqVWrljLGykbasmWL7N69W93va8oZMHBwuiIt6tNImpz4neAeYCvhzA0+wEzRQ8APAY5erwpuUS4J94EDP8muXTulbNlyUmT\/dtOoYnvTILBq1Sp58MEHZfCMF+K7BQ0BXrx4cSCDwTQBTmYJfdnxn8i2kb9uRSNRgAMZOhlVkosCnBFQhj5Efg0lTrPZsbeCxrV\/ixYtEkS8CiJhq7tOnTpBVKVVRzqCkxlijau0ircgaaEbfiZO0OFzkM0WkN9soht+2enm5\/BbWHAL0lpBR70D2W5fOoKdAnx+rZPkiXKL8l2pl9dppORdMSLbTWX5GSDACToD0Ax6hPwaRFYGTU03P2dQZKCPUIDTwK1DsNMQa13p0flcWngLUqBj2lVlnKBdwWVcZvJrHGWuGqwzP7sqMODMFGAfBNhpiPXq4Tul3LZ3EiXDSbxE\/eYBU8vqdBDgBK2Dkrl5yK+53Om0nAKsg5LBeXQIdhpivbe5Zb4e11xYYKwTg9Exv+mcoM3nsKAekN9486szP0cZAa6AfVgB28+BKx3eLou3dU6USgvoKA9\/EU7Q0ebHa+vIr1cEo\/08BTja\/HhunQ7BdgE+5+e1MmNH\/0S9xes3l8ojl3luBwvIDgKcoLODa1RKJb9RYSI77dCZn7NTsz+lcgXswwoYRViGWL2+f1x67Xs8UWqp5t0FRlhM0USAE3Q0efGrVeTXLySjWQ4FOJq8+NYqXYItQ6wRe8ZL2x\/\/kaifLki+UZGVgjhBZwXWyBRKfiNDRVYaojs\/Z6VyHwrlCtinFbBliDVjx+1yzs8fJEqlBbQPozSLRXCCziK4ESia\/EaAhCw2gQKcRXCjULQuwdY58KJtnaXy4aMxZxEDOtfuX40Cb7pt4ASti5SZ+civmbzptlp3ftYtL+h8XAH7tAJGMfVveZoW0EGPYI\/1cYL2CGDEHye\/ESfIY\/MowB4BjPrjbggeNXmudH27R6JLtICOOrt0Q4o+Q95aSAH2hl\/Un3YzP0exL1wB+7gC3rtgVL4Y0BTgKA75\/G3iBB19jry0kPx6QS\/6z1KAo8+Rpxa6IdgpwDNLd5M+E6bLKWWyc4eyp47xYYUAJ+h4DwTyG29+3czPUUSCK2AfV8BbR7bIdwlD7\/KTpWf3q+XqxhWjyD3bRAGO\/RigAMebYgpwvPkVNwR\/3aeGHN65KYFI20rzpfkfz5KpV9eLOUrmdo8TtLnc6bSc\/OqgZG4eN\/NzFHvJFbBPK2AILwTYSluLVJR2learf64d\/kduQ0dx9HMFHFFW\/GsWBdg\/LKNYEgU4iqz42CZdgg98tFy2jWyRqPn9Yg2ld\/kH1L+xBc1VsI+k+FgUJ2gfwYxgUeQ3gqT42CTd+dnHKn0tiitgn1bATgOsxSe0kVFlBidK3zP5qDj7yiAL84QAJ2hP8EX+YfIbeYo8NZAC7Am+6D+sS7DTAAsW0DNP7Jbo4OA21WVwm6Nb1NHveW60kBN0vHkmv\/HmV3d+jioKXAH7tAJ2CvDyS56Wgf+pkK90ngVH7zXgBB09TvxsEfn1E83olUUBjh4nvrZIh2CnARYaUHPhEbFuSLIa1OmcCjKjy+99bR8L84YAJ2hv+EX9afIbdYa8tU9nfvZWQ3af5grYhxWwU4Bx+QIuYbAuaLBXsahPI2la+6TsssrStRHgBK0NlZEZya+RtGk3mgKsDZWZGXUIdlpA20NQtp26Rv71xXf5Ok8Rjs5Y4AQdHS6y0RLymw1Uo1OmzvwcndYe2xKugH1YATstoPM6jZS8K0aokld+\/p20m7YmXy0ITYnzYKbwEeAEHT4H2WwB+c0muuGXTQEOn4OstkCHYKcBVqWRy6RE\/eaJdqUS4SlX1eN2dFbZS184J+j0GJmcg\/yazF76tuvMz+lLCS8HV8AeV8DJDLBw\/otzYHsasfgLeXjZ18fURvek8AY\/auYEHS7+2a6d\/GYb4XDLpwCHi3\/Wa09HcCoDrGQNu\/qx\/8g\/Ptp9zK+wJQ0fYV7akHU6j6mAE3TwmAdZI\/kNEu3g60o3PwffInc1cgXscQVckAFWsqLve3WjTHz96IUN9jwQ4r7Nq0m9iiW5Ne1uHGecmxN0xtAZ8SD5NYKmjBtJAc4YOjMeTEdwQQZYqXqYzD0pWV7rHuHWvz9ZGlUrnchSLe\/X+4Xt9wzzzuHMxhMn6MxwM+Up8msKU5m1M938nFmpwT3FFbDHFXA6A6xUxX+95yd5+t1tMv4fyVfD2RwCFOv86B4+fFiKFCmSTchZdogIkN8QwQ+g6i1btsjL3atJkyZNAqjN3yoowB4EOFUELDcUQYhXbNgr\/Z79xM1jzEsEiAARIAK\/IbDkihIU4DiOhoK2OJznv1YErExxgBhXzSsml05be0zwjkzL5HNEgAgQgbgjQAGOKcMFCbDz\/NceAcsvOOBDvHnvT6q4r\/ccEIg0fpCs\/\/\/1d7\/+HxMRIAJEINcQoADHlPGCBNh5\/luu7xwp1bx75JGgWB+lCEY6O3fulHLlyknx4r8atzHFBwHyGx8uk\/Vk1ap3ZNCgO+SZWQ9xCzqOVBckwF92KpSvy84IWHHEI259opVs3BjN3x\/yG29+aQUdb34lFcFuAnDEHCKju8cJ2mj60jae\/KaFyOgMFGCj6Uvf+FQEB3H+m751zOEVAU7QXhGM9vPkN9r8eG0dBdgrghF\/PhXBzvNf+w1IEe8Sm2dDgBN0vIcD+Y03vxTgePObcgua57\/xIJ4TdDx4TNUL8htvfinA8eY3qQDz\/Dc+pHOCjg+XyXpCfuPNLwU43vwmFWCe\/8aHdE7Q8eGSAhxvLpP1jgIcc86TEWyq\/2\/MqcqoexTgjGAz5iHyawxVGTWUApwRbOY85CQ4Wfxn+v+aw6ezpZygzeVOp+XkVwclc\/NQgM3lTqvlToL9jv+s1QhmyhoCnKCzBm0kCia\/kaAha42gAGcNWncFHzlyRD777DN5\/PHH5fXXX5fdu3dLqVKlpFmzZtK1a1dp3LixFC5c2F2hIsecAdP9yDWEkX6AE3Sk6fHcOPLrGcJIF0ABjgA9uO9z\/vz5MmbMGMHfk6Xrr79ebr\/9dilRooSrFtsJPrdWRfm6T418z3P72RWckcvMCTpylPjaIPLrK5yRK4wCHAFKXnjhBRkyZIgS36ZNm0qfPn2kTp06snHjRpkxY4YsXbpUtfLWW2+Vvn37urp83U5ww1I\/ybaRLRI99nr9YASgy\/kmcIKO9xAgv\/HmlwIcMr\/ffPON9OrVS9avXy\/t27eXUaNGScmSJROt+vnnn2XKlCkydepUdePN3LlzpV69etqtthNc6\/0Zsn\/53MSzuPkINyAxmYsAJ2hzudNpOfnVQcncPBTgkLl79tlnZejQoUpcH3vsMWnQoMExLdq+fbv07t1bPvzwQ7UNffPNN0uhQvlvMkrVDYvgJx4cK5WfuIrbzyHz7Xf1nKD9RjRa5ZHfaPHhd2sowH4j6qK8AwcOqK3nxYsXS6tWrWTy5Mn5Vr9WUTDQmjBhgkyfPl3OOussmTlzppQtW1arJovg2cN7yymvDuH2sxZq5mTiBG0OV5m0lPxmgpo5z1CAQ+Rq27Ztct1118knn3yiznYHDBiQsjWvvPKK9OvXT04++WSZM2dO0pVysoctgifV2SY4A7YSt59DJN7HqjlB+whmBIsivxEkxccmUYB9BNNtUXA76tGjh0CIx44dK1deeWXKIt59913p0qWLMtSaNWuWtGzZUqs6EHx79yvkqTM258tP62ct+CKfiRN05Cny1EDy6wm+yD9MAQ6RIgt8NAFuSE2aNEnZmg0bNiix3rp1a1qxtheCOt6640\/SqcL3if8uXr+5VB65LMSes2q\/ENi0aZPyHe\/WrZtUr17dr2JZTkQQIL8RISJLzaAAZwlYnWKDEuDyk\/7I1a8OIQbmMf0FNhDyQJtMfgOFO\/DKTOe30BFYKBmawhDg7QeLSOGhKwxFjM12IgA3tkGDBikf8fPOO48AxQwB8hszQh3dsfhNtwMaVRQowGmY2bJli+yc2kNO\/HypQHwnbConH+wvHlU+2S4iQASIQE4hgA9neLlUrVrVuH5TgDUogwjjh4kIEAEiQASihQCE10TxBYoU4GiNJbaGCBABIkAEcgQBowXY7oY0adIkFYoyVbKfFz\/xxBMqZjQTESACRIAIEIGwEDBagHfu3Cm45QghJvv3769CTKZKL774ogrU4TYQR1jEsF4iQASIABGINwJGCzAuWrjzzjsF4tq2bVsZN25c0usGvYSijDf97B0RIAJEgAiEhYDRAgzQdC5jsN+YhNCVgwcPdnUlYVjksF4iQASIABGILwLGC7BdXC+77DIZOXKklC5dOsGY8zrCVDcmxZdi9owIEAEiQASiiIDxAozt5eeee06GDRum4jw3btxY+vTpI7\/\/\/e\/l66+\/lhkzZsjSpUsV9gi2gEsbihQpEkUu2CYiQASIABHIIQSMF2BwBeHFDUdwxsbfkyUYa+Eu4BIlSuQQvewqESACRIAIRBWBWAgwwMVKGG5JCKz\/+uuvy+7du6VUqVLSrFkz6dq1q1oZFy5cOKo8sF1EgAgQASKQYwjERoBzjDd2lwgQASJABAxHgAJsOIFsPhEgAkSACJiJAAXYTN7YaiJABIgAETAcAQqw4QSy+USACBABImAmAhTgJLz9+OOPsmTJEnn66adl3bp1yrK6bt26Aj\/jjh07SpkyZcxkOwdbbXdTu+eee+TKK69MiwJ8x1esWCGIGb5q1So5dOiQVKlSRUVb69y5s7E3r6TtuCEZfvnlF3n33XflySefFMR4h8ElXAsbNGigOLr00kvTvqPkOLpkp+IGBrXdunWTOnXqSKFChQrsgCn8UoAdNG7evFnFlX7\/\/feTElyjRg2ZOHGiNGrUKLojmC1LILBy5Url+71\/\/34ZO3ZsWgHeu3ev3H333fLKK68kRRGW9XB3u\/jii9NOAqTBfwQQ\/33UqFHqAzlVqly5suByFtwTmyyRY\/958atEBFYaOnSo4L1NlvChBZfSfv36pXQpNYlfCrCNZbzcCNaBr2pMtDfeeKO0a9dO5XjjjTfk4YcfVl\/bEOGZM2dKrVq1\/Bp3LMdnBLDyxUuMjylwhpROgA8cOCB33XWXvPDCC2pFha9trHhPOOEEtRKePHmyfPXVV2psIMBLkyZNfG41iysIgR9++EFGjBihYr8jYUeqS5cu6n3ct2+fvPzyy8oN0XJBnDt37jEfyuQ4umMMwon5F++t9f516NBBKlSoIBs3btQKqmQavxTg38YjJmxsOeLrGhPs1KlTj7mycP369erGJQwGxpSO7ouMI4TZs2erDyZ7YJZ0Avzmm2+qjy6ke++9Vx032Le68HWOG7VWr14trVu3VjshJUuWjC4QMWsZVr233Xab4hQrIPw4o9qtWbNG\/f\/WrVuTckSOozsorLj+4DTZ+2f\/ACtXrpzgA6tevXr5OmQavxTg3+jbtWuX9OrVS9auXSvdu3dXtywlC1kJkUa8aXx1I6509erVozuic6xlmJj\/+c9\/Cs56sVJ1poIEGF\/OQ4YMkcWLF0urVq3UFiY+xJzprbfekhtuuEGKFSumJoBzzjknx1AOp7t2frDzgA\/kvLy8YxqDD+kpU6bIAw88cMzVo+Q4HO50aoW4YrcKYYMLev9gk9OjRw+1y+G8A95EfinAv40ObDtfe+216ut61qxZ0rJly6TjZsOGDWoA4As73YpKZ+Axj38IgENsGVsJ0c+uueYaGT9+fFq+7LwOHz5cevbsmbRh9g81rJYHDRrEs2D\/KExZEmwz8N59+eWXae\/+to+D+fPnJ44KyHEARGVYBY7\/sLjB9nObNm0SO1HO4uwcwrYDO1JWMpFfCvBv7GHLcsyYMVKtWjW1ssEKN1nCWRO2oTFQEOISZ4a83CHDt87nx6yJ97jjjlPbkNaHks4HE7auYNwBLp966ikVujRZOnjwoDLSWrBggfpSx7kwt6F9JjJJcZhccf4LIcb7V5A1eyoBJsfZ5ynbNRS0AjaRXwrwb3GksZrFyrdhw4ZqazmVqxHM27E9DUMQTsDZft3clQ\/L9eXLl6udDJwRIenuWDzzzDOKV1jQ4mKP0047LWXl2PrCFihc0zBWKlWq5K6hzJ1VBKyPaRjPYVVleSyQ46zCnvXCMffed999Mm\/ePPV+O8+ATeSXAiwi9rMD+BGOGzeuwFuTsKUJK9h0Yp31EckK0iKgK8BuOLWMRXTEOm0DmcFXBLZv3652P\/Ax5jwrJse+Qh1YYdh1+ve\/\/y2PPvqoLFu2TO1SJTPSMpFfCjAFOLAXKYyKKMBhoB5OnbDfgOU7fpCcNhomTtDhIBmNWvfs2aO8TT744INEg0499VR1BHThhRcec7udifxSgCnA0XjbstQKCnCWgI1YsfaIZxDi9u3bK5dC+\/m8iRN0xGAOtDl2wzt7xdjZGDx4sJx55pn5DCBN5JcCTAEO9KUKujIKcNCIB18fxBc2GYighLChf\/jDH+Shhx6S8uXL52uMiRN08GhGp0YYvG7ZskUqVqyYNNCKM1aDifxSgCnA0XnjstASCnAWQI1QkYgLvXDhQuWNgJUv7DJgJFezZs1jWmniBB0hqCPRFHswJAThQERCxGlHMpFfCjAFOBIvVrYaQQHOFrLhlwur2EceeUSmTZumxBdBUeAWBlfCZMnECTp8lKPXAisYElqG8\/5LLrmEAhw9mvRbhC0syw0JLzG+qpJF2UGJdotpEA+Labg7MEUTAV0BtlwY0vmB27+0zzjjDOW6Zrk8RROBeLYK25OIePb888+rDsIlEJaxBXFBjuMxFj777DPl479t2zZ10YoVjMNEfrkC\/m1MWuRh6wp+oKm+ohEwHCEr4ebAQBzRf6F1BRiBVeA\/7PQddfaQfuDhc75jxw4ZOHBg4sacyy+\/XG1Bly5dusDGkePwufOjBfZ3unfv3sogC8lEfinAv40I3C+Km1WQsMWR6qYb3Qndj4HGMrwjoMsXLthADHBYXhYUYpShKL1z4qUEiO8tt9yiLsRAwgoIkbEQmztdIsfpEArv959++qmKsY8Qv0vg2ZoAAA2oSURBVDgqKOimMVy4gY9lXLpiDxtrIr8U4N\/GnM5lDPYbk1LdxhHeEGbNyRDQFWDdQO7221YQGAD+iEzBIGC\/DQfBGBCHG1uRuqFgyXEwPGVSi\/09LegyHJz1P\/jggyoSnXO3ykR+KcC\/jRbndYRwY7jgggvy+ZnZLfAKGiSZDEA+kx0EdAUYtdvFFeeLnTp1yufsb7+OsKAbW7LTk9wuFe\/n9OnTZcKECQqIVNcRpkOJHKdDKJzf42hn9OjR8vTTT6tbyJLNvxgDr7\/+uvrw2r9\/v+CuYLynJUqUSDTaNH4pwLbxZt\/eQkB\/nPXifAl\/f+ONN5TFHa7BwkUNMNSqVatWOKOVtWoj4EaAnSssCDBuUzr55JNl1apVysIW1xxigkAo0oK2ybQbyIxaCNh5ROxtTMJly5ZN+ywiJ9ntOchxWshCywDjqptuukndt445F7sbl156qVSoUEG+\/fZbefLJJ5XLGVbBtWvXVtdO1qlTJ197TeOXAuwYbrjuDFZ19vBn9iwQX1zEbgV4D220smItBNwIMArEtWjDhg1T95ImSxBfrMIuvvhiXkOoxYA\/mawLFtyWluw8nxy7RTG4\/BBhGFWlmn\/REnz44pw4laGsSfxSgJOMLRzuL1myRG2H4PorfHHh5pvLLrtMOnbsmPKmpOCGKWvSRcCtAKNcbIetWLFCGeNh5YvoSnD2x0UduG+4atWqutUznw8I2K+AdFtcKoM6cuwWyeDy4ywXH8BY7VrvHz58Iby4hvL8889Pa3RnCr8U4ODGFWsiAkSACBABIpBAgALMwUAEiAARIAJEIAQEKMAhgM4qiQARIAJEgAhQgDkGiAARIAJEgAiEgAAFOATQWSURIAJEgAgQAQowxwARIAJEgAgQgRAQoACHADqrJAJEgAgQASJAAeYYIAJEgAgQASIQAgIU4BBAZ5VEgAgQASJABCjAHANEgAgQASJABEJAgAIcAuiskggQASJABIgABZhjgAgQASJABIhACAhQgEMAnVUSASJABIgAEaAAcwwQASJABIgAEQgBAQpwCKCzSiJABIgAESACFGCOASJABIgAESACISBAAQ4BdFZpFgK\/\/PKLfPzxx5KXlydVqlQxq\/GGt\/abb76RvXv3yu9\/\/3spXLiw4b1h84lAfgQowBwRsUFg\/PjxMmPGDM\/9adWqlUyePFlKliwpR44ckenTp8uECROkSJEi8sQTT0iTJk081xHFAuz44UOjbNmyqpnFixeXUaNGSZ06dQJt9po1a6R79+6yf\/9+ufHGG2XQoEFSqFAh39qA8kePHp0ob8uWLbJ792717+HDh0vPnj19q4sFEYFkCFCAOS5igcCBAwdkyJAhsnjxYqlUqZK0a9dOGjduLEWLFlX9w0pq4sSJaoI9+eSTZeDAgYnV7MGDB+Xdd9+VRYsWybZt28QuwD\/\/\/LPceeed8uKLL6pyRo4cKddee20sMHN2whLg1q1by9ixY9WKP8wEzAcMGKCagDaBP3wUZSthjPTq1UvWr19PAc4WyCw3HwIUYA6IWCBgCfDWrVvV6rVatWr5+rVhwwbp0aOH4PcQ4Dlz5kiDBg3y5cFW59ChQ+Xrr7+Wxx57TAk5EsT5jjvukNNOO01GjBgR221oS4C7du0qd911l1rxh5l27typVt4fffSR+vC54IILfF0BO\/tm\/4jjCjhM5nOnbgpw7nAd655ak+dFF10kbdu2PWai1hFgAPT+++\/LmDFjZNq0aQkBjjVwts5FTYCDxp0CHDTirI8CzDEQCwQwed57773SrVs3tVJ1Jl0B3rVrl1p1YdvZWgHHAiCNTlCAjx5jcAWsMWCYxTMCFGDPELKAKCDw448\/yuuvvy4tWrSQE088MWMBxnnw0qVLlaFVmTJlotC1wNpAAaYABzbYWJFCgALMgZATCOiugHMCjBSdpABTgHN5\/IfRdwpwGKizzsAR8CLAP\/zwg\/zf\/\/2fMtyClSxW2fYEP2FY0L788suyceNGufvuu+W4446TFStWyLx582Tt2rXy008\/KUteWGpbvsRYtf\/jH\/+Qv\/3tb7Jq1SopXbq0tG\/fXm6++Wb191TJeu6ll15KlH3WWWfJddddJxdeeKEUK1YsI3x1Bfjw4cMCPGFxDjctWCrDWnzJkiXy3HPPqTbBdenyyy+Xfv36Jayp9+zZI2gznlu3bp1UqFBBGcZ16dIlaZtRz2effaZcv4AHDOGchmFwUfrnP\/+prNTRjnr16ikjOriOYUfk+++\/l2bNmimr5ho1ahSIC8+AMxo2fMgDAhRgD+DxUXMQcCvAEJb33ntPnn32WSUYhw4dUp2dNWuWtGzZUv0d29VPP\/20EtAPPvhA\/R9cmAYPHiwQM\/iVYisbrk9\/\/\/vfVRkQypkzZyrBhlDXqlVLCTLqggEYRAdGZOPGjZMSJUocAzBcZPr37y9NmzaVDh06KItuiDfqg4X3X\/7yF3UWnmwbPh1b6QQYwv\/II4\/Ia6+9Jp9\/\/rkqDhbTV111lQwbNkz9+7zzzlOiiY8Pqy\/33XefvP3226pPwKNUqVKycuVKZd2M5PTxRT8WLFggCxcuVG5hVj12y+zVq1fLo48+Km+99ZaqBziAGwj7\/fffL1WrVlUfPfggQqpcubL6gEpmH2DhQgFON0L4e78RoAD7jSjLiyQCmQgwVldY2cElCaLpFOBkEzdWvv\/v\/\/0\/tYpt1KiRit4EMcfKECIFsYBIQlxvv\/12JQxIyIOVHgzAIFD4e8OGDfNhCWG76aablB8yVo321eCbb76phAzlI3gFjMjcuhGlE2CrMTBUw04AVronnHCCtGnTRnr37i21a9dW1udow8MPP6x+sBrHRwkM2uwre+SBMM+dO1dq1qyp\/oRoIkEI8fPvf\/9bYQThT+YaZccM9UDcIbA33HCDlCtXTmG6fPlyueWWW1QZaCM+jlIlCnAkX91YN4oCHGt62TkLAbcCbEfunXfekc6dO6cUYPzCEi9sc2Il5tzuhE\/r9ddfLx9++KFavU6aNEmJhD1BYLEli1WfM+AHPgQQtQmih48BZ0hMuyimEvB0o0FXgO1CpdMXbEWj7c4VPVbBVlATrE6xfW5P2LLGtjp2F1L5Jtu5wQoZPNm34CG82PZ\/5ZVXVPkPPfSQ+sBJlijA6UYIf+83AhRgvxFleZFEICgBrlu3br4gHslWyfZIW6kEx7laQ9hEiBXOiLF17VzdYjsc\/4+tW6RMInZlIsA6fUklnnZOEHnryiuv9CTA9uMBe0FWv1Jxk4wjuiFF8jWOXaMowLGjlB1KhoBpAuwUrdmzZ6sAIbopk2hWYQpwMsFzuwKmAOuODuaLCgIU4KgwwXZkFQGTBRjnpffcc4+yqIYBFs5Ss5EowHRDysa4YpmpEaAAc3TkBAImC7D9bDKTla0uwRRgCrDuWGE+fxCgAPuDI0uJOAJxEeBUZ65+wE8BpgD7MY5Yhj4CFGB9rJjTYARMFmC408BICWeccOeBxXC6u3lhdQ1r4IICejjppABTgA1+xY1sOgXYSNrYaLcImCzA6CsCguCqRCS45sCfNZWf744dO9SVjLfddptUrFhRGyoKMAVYe7Awoy8IUIB9gZGFRB0BuwCjrfPnz1eBG3SSGz9gv9yQnGe9mzZtUsKLyE4QXkS76tix4zHXLu7bt08F8zj77LOVW4+bYBwUYAqwzvvAPP4hQAH2D0uWFGEE7EEu0MxULivJuhAFAbZHfbLaiPNgRHlCNC0YaiH844wZM6RBgwYq8EVeXp4rRijAFGBXA4aZPSNAAfYMIQuIOgIQLwSosLZw0d5rrrkmaUCLZH1BWEgEtkAaNGiQCgdpT3YrZUTAQujK6tWr58uzd+9e6du3r0DMnaEXrYz2Vfoll1yiYicj1KOVEA0LsZinTZumwj0mS1dffbUMHDjQtfiiLF0Btkfdwi4C2uSMPW1fsSNa1pQpU445jwYWCC6CviBKGPhBKEsrIe41wmriPDtVGXZuENoScantCTsCcN1CqM505+eMhBX1Nzl+7aMAx49T9ug3BLZv3y7PPPOMvPHGG4nA\/3ZwEM4R4Ql79uypRNGeIHa4dACXCDz\/\/POJyxiwpYvYx7gwAbfs4O5g5MOtRpYonnrqqdKpUyf1U7RoUXnhhRdk2bJl6gICKyEPxAWRrfDcU089pW4T+vTTT1UWqx7coIQfK7wiPiZwiQFu+8GKF\/GqcREBbmjC5Qznnnuuq21ne5\/TCTDEDBdToK\/2viBmNcJNtmvXTsVcBl7A3LqgAnUgD7bV8WGBM2r4NCM8pHXZAmJoowxge\/zxx6sbjnC71FdffZVoIsoA7siHtkB8cREGLrtAQhl\/+tOfVCxoHAVAdJ988sljcEc78KHiDI1JAebUETQCFOCgEWd9RCCiCKQT4Ig227dmUYB9g5IFaSJAAdYEitmIQNwRoADzDDjuYzxq\/aMAR40RtocIhIQABZgCHNLQy9lqKcA5Sz07TgTyI2AJMO4zxl2+bq2oTcfTbijH25BMZ9OM9lOAzeCJrSQCWUfAEmBUBMOuqlWrqjphrATfYlx2H6eEKx7hrmUlWHd\/88036p8U4DgxHd2+UICjyw1bRgSIABEgAjFGgAIcY3LZNSJABIgAEYguAhTg6HLDlhEBIkAEiECMEaAAx5hcdo0IEAEiQASii8D\/B6fnohBwnj2xAAAAAElFTkSuQmCC","height":193,"width":320}}
%---
