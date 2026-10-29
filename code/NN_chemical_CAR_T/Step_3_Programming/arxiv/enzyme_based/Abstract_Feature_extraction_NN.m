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
% brendaKineticsTable('Selected_Full_Chain_1.txt','kinetics_BRENDA.xlsx', ... %[output:group:0077ef08] %[output:66d59597]
%         'BrendaFlatFile', which('brenda_2026_1.txt'), 'KcatSource','BRENDA else SAxMW'); %[output:group:0077ef08] %[output:66d59597]
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
% K_EV_AC = 10*vesicles(row,1);            % activation constant in DATASET EV units (assumed) %[output:812f1865]
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
m_cAMP = cAMP_ref*ones(N,1);
%m_AMP = AMP_ref*ones(N,1);
% m_ATP = ATP_ref*ones(N,1);
% m_glu = glu_ref*ones(N,1);   m_pSer = pSer_ref*ones(N,1);
m_in = m_in_target*ones(N,1);
m_Ff_1  = m_Ff_ref*ones(N,1);    m_Sf_1  = m_Sf_ref*ones(N,1); m_Sd_1  = m_Sd_ref*ones(N,1);  m_Fd_1 = m_Fd_ref*ones(N,1);  m_d1 = m_Fd_ref*ones(N,1);
m_F2 = m_F2_ref*ones(N,1);   m_S2 = m_S2_ref*ones(N,1); m_d2 = p2_ref*ones(N,1);  m_N2 = p2_ref*ones(N,1);
m_c = m_c_ref*ones(N,1);     m_I = m_I_ref*ones(N,1);   m_out = m_out_ref*ones(N,1);

for n = 1:N-1
    dt = t(n+1)-t(n);
    % Block 0
    v_EV = V_EV*EV(n)/(K_EV_AC+EV(n));                               % EV count -> cAMP (assumed)
    v_a1 = V_a1*m_cAMP(n)/(Km_a1_cAMP+m_cAMP(n));                    % cAMP         -> AMP
    % v_a2 = V_a2*m_AMP(n)/(Km_a2_AMP+m_AMP(n))*a2_fraction;           % AMP          -> ATP
    % v_a3 = V_a3*m_ATP(n)/(Km_a3_ATP+m_ATP(n))*a3_fraction;           % ATP          -> L-glutamate
    % v_a4 = V_a4*m_glu(n)/(Km_a4_glu+m_glu(n))*a4_fraction;           % L-glutamate  -> 3-P-L-serine
    % v_a5 = V_a5*m_pSer(n)/(Km_a5_pSer+m_pSer(n));                    % 3-P-L-serine -> L-serine

    m_cAMP(n+1) = m_cAMP(n)+dt*(v_EV-v_a1);        
    m_in(n+1)   = m_in(n)  +dt*(v_a1-(k_Sf_1+k_Ff_1)*m_in(n));
    % m_ATP(n+1)  = m_ATP(n) +dt*(v_a2-v_a3);
    % m_glu(n+1)  = m_glu(n) +dt*(v_a3-v_a4);
    % m_pSer(n+1) = m_pSer(n)+dt*(v_a4-v_a5);
    % m_in(n+1)   = m_in(n)  +dt*(v_a5-(k_Sf_1+k_Ff_1)*m_in(n));

    k_I_1=5;
    % ---- Differentiator 1 ----
    m_Ff_1(n+1) = m_Ff_1(n)+dt*(k_Ff_1*m_in(n)-k_Fd_1*m_Ff_1(n));                    % fast delay pool
    m_Sf_1(n+1) = m_Sf_1(n)+dt*(k_Sf_1*m_in(n)-k_Sd_1*m_Sf_1(n));                    % slow delay pool
    m_Fd_1(n+1) = m_Fd_1(n)+dt*(k_Fd_1*m_Ff_1(n) - k_Fr_1*m_Fd_1(n) - k_I_1*m_Fd_1(n)*m_Sd_1(n));
    m_Sd_1(n+1) = m_Sd_1(n)+dt*(k_Sd_1*m_Sf_1(n) - k_Fr_1*m_Sd_1(n) - k_I_1*m_Fd_1(n)*m_Sd_1(n));
    m_d1(n+1)   = m_Sd_1(n+1);

    % m_Ff_1(n+1) = m_Ff_1(n)+dt*(k_Ff_1*m_in(n)-k_Fd_1*m_Ff_1(n));                    % fast delay pool
    % m_Sf_1(n+1) = m_Sf_1(n)+dt*(k_Sf_1*m_in(n)-k_Sd_1*m_Sf_1(n));                    % slow delay pool
    % m_Fd_1(n+1) = m_Fd_1(n)+dt*(k_Fd_1*m_Ff_1(n) - k_Fr_1*m_Fd_1(n)         - k_I_1*m_Fd_1(n)*m_Sd_1(n));
    % m_Sd_1(n+1) = m_Sd_1(n)+dt*(k_Sd_1*m_Sf_1(n) - (k_Ff_2+k_Sf_2)*m_Sd_1(n) - k_I_1*m_Fd_1(n)*m_Sd_1(n));
    % m_d1(n+1)   = m_Sd_1(n+1);

    % P = m_Fd_1(n)+dt*(k_Fd_1*m_Ff_1(n)-k_Fr_1*m_Fd_1(n));                         % fast rail
    % Q = m_d1(n)+dt*(k_Sd_1*m_Sf_1(n)-(k_Ff_2+k_Sf_2)*m_d1(n));                % slow rail; removal = load of D2
    % D = P-Q;                                                           % canceler, implicit over dt
    % qb = 1-dt*k_I_1*D;
    % qd = sqrt(qb^2+4*dt*k_I_1*P);
    % if qb > 0
    %     P = 2*P/(qb+qd);
    % else
    %     P = (qd-qb)/(2*dt*k_I_1);
    % end
    % m_Fd_1(n+1) = max(P,0);
    % m_d1(n+1) = max(P-D,0);
   

    % ---- Differentiator 2 (fed by the m_d1 it consumes) ----
    m_F2(n+1) = m_F2(n)+dt*(k_Ff_2*m_d1(n)-k_Fd_2*m_F2(n));             % fast delay pool
    m_S2(n+1) = m_S2(n)+dt*(k_Sf_2*m_d1(n)-k_Sd_2*m_S2(n));             % slow delay pool
    kI2_now = k_I2c*(m_I(n)/(KM_Cr+m_I(n)))/satI_base;                 % r_I2 at the present m_I
    P = m_d2(n)+dt*(k_Fd_2*m_F2(n)-(k_c+kI2_now)*m_d2(n));             % fast rail = output; removal = counter
    Q = m_N2(n)+dt*(k_Sd_2*m_S2(n)-k_Sr_2*m_N2(n));                    % slow rail
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

    pools = [m_out(n+1),m_in(n+1),m_Ff_1(n+1),m_Sf_1(n+1),m_Fd_1(n+1),m_d1(n+1),m_F2(n+1),m_S2(n+1),m_d2(n+1),m_N2(n+1),m_c(n+1),m_I(n+1)];
    assert(all(isfinite(pools)) && all(pools >= 0),'Invalid step at n = %d; reduce dt0.',n);
end

%[text] Diagnostics
% The identity behind the design, computed by plain filtering of m_in:
% m_P1 - m_d1 = w_1*(m_in delayed by tau_F1 - m_in delayed by tau_S1), smoothed with tau_r1.
x_F = m_in;  x_S = m_in;  ideal_1 = zeros(N,1);
for n = 1:N-1
    x_F(n+1) = x_F(n)+dt0*(m_in(n)-x_F(n))/tau_F_1;
    x_S(n+1) = x_S(n)+dt0*(m_in(n)-x_S(n))/tau_S_1;
    ideal_1(n+1) = ideal_1(n)+dt0*(w_1*(x_F(n)-x_S(n))-ideal_1(n))/tau_r_1;
end
diff_1 = m_Fd_1-m_d1;
fprintf('\nD1: m_P1 - m_d1 against the ideal subtraction: largest gap %.3g mM (signal up to %.3g mM)\n', ...
    max(abs(diff_1-ideal_1)),max(abs(ideal_1)));
fprintf('D1 output m_d1: rest %.3g, peak %.3g at %.3g min, last %.3g mM;  m_P1 peak %.3g mM\n', ...
    m_d1(1),max(m_d1),t(find(m_d1 == max(m_d1),1))/60,m_d1(end),max(m_Fd_1));
fprintf('D2 output m_d2: rest %.3g, peak %.3g at %.3g min, last %.3g mM;  m_N2 peak %.3g mM\n', ...
    m_d2(1),max(m_d2),t(find(m_d2 == max(m_d2),1))/60,m_d2(end),max(m_N2));
fprintf('Counter m_c: rest %.3g, peak %.3g (%.3g times the rest) at %.3g min, last %.3g mM;  m_I lowest %.6g mM\n', ...
    m_c(1),max(m_c),max(m_c)/m_c(1),t(find(m_c == max(m_c),1))/60,m_c(end),min(m_I));
target_N = w_N*m_c+b_N;
fprintf('Neuron m_out: rest %.3g, peak %.3g at %.3g min, last %.3g mM;  largest gap to w_N*m_c + b_N: %.3g mM\n', ...
    m_out(1),max(m_out),t(find(m_out == max(m_out),1))/60,m_out(end),max(abs(m_out-target_N)));
fprintf('Largest pool / K_M (first-order kinetics need << 1): m_F %.2g, m_S %.2g, m_F2 %.2g, m_S2 %.2g\n', ...
    max(m_Ff_1)/KM_Fd_1,max(m_Sf_1)/KM_Sd_1,max(m_F2)/KM_Fd_2,max(m_S2)/KM_Sd_2);


%[text] Figures
% 1. Input.
figure;
plot(t/60,m_in,'LineWidth',2); grid on;
xlabel('Time [min]'); ylabel('L-serine: m_{in} [mM]'); set(gca,'FontSize',fontsize+4);

% 2. Differentiator 1: the two delay pools.
figure;
plot(t/60,[m_Ff_1/m_Ff_1(1),m_Sf_1/m_Sf_1(1),m_in/m_in(1)],'LineWidth',2); grid on;
xlabel('Time [min]'); ylabel('relative to rest'); xlim([0,min(30,t(end)/60)]);
legend('m_F (fast delay pool)','m_S (slow delay pool)','m_{in}','Location','best'); set(gca,'FontSize',fontsize);

% 3. Differentiator 1: the two rails, and the falling slope of the input.
figure;
yyaxis left;  plot(t/60,[m_d1,m_Fd_1],'LineWidth',2); ylabel('Concentration [mM]');
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
% 8. Enzyme Concentrations Bar Plot
figure; %[output:49e91d33]
% Combine the calculated enzyme values in uM

% Create a cell array matching the labels for all enzymes
enzyme_labels = {
    'r_F'; 'r_Fr'; 'r_S'; 'r_Sr'; 'r_I'; 'r_Fx (new)'; ...
    'r_F_2'; 'r_Fr_2'; 'r_S_2'; 'r_Sr_2'; 'r_I_2'; 'r_Sx_2 (new)'; 'r_c'; 'r_I2'; ...
    'r_Nw (Neuron)'; 'r_Nb (Bias)'; 'r_Nr (Removal)'; 'r_I1 (Counter)'; 'Block 0 (Total)'
};

% Create the bar plot
bar(E_all*10^3, 'FaceColor', [0.2 0.6 0.8]); %[output:49e91d33]
grid on; %[output:49e91d33]
set(gca, 'XTick', 1:length(enzyme_labels), 'XTickLabel', enzyme_labels, 'TickLabelInterpreter', 'none'); %[output:49e91d33]
xtickangle(45); % Rotate labels so they don't overlap %[output:49e91d33]

% Formatting labels and titles
xlabel('Enzyme / Reaction Link'); %[output:49e91d33]
ylabel('Enzyme Concentration [\muM]'); %[output:49e91d33]
title('Calculated System Enzyme Concentrations (\muM)'); %[output:49e91d33]

% Dynamically adjust font size if defined in parameters
if exist('fontsize', 'var') %[output:group:1c204f29]
    set(gca, 'FontSize', fontsize); %[output:49e91d33]
end %[output:group:1c204f29]
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
%[output:66d59597]
%   data: {"dataType":"text","outputData":{"text":"Parsed 27 reaction occurrences across 1 motif(s).\nNo BrendaFlatFile given; using the local BRENDA flat file on the path:\n  D:\\code\\bioelectrodynamics\\code\\database\\brenda_2026_1.txt\\brenda_2026_1.txt\n(SOAP is bypassed.)\nLoaded BRENDA index cache: 18 EC(s) from D:\\code\\bioelectrodynamics\\code\\NN_chemical_CAR_T\\Step_3_Programming\\enzyme_based\\brenda_ec_index.mat\nIndexing local BRENDA for 2 new EC(s) (18 already cached)...\n  Done: 2 new EC(s) indexed (20 total cached).\nSaved BRENDA index cache (20 EC(s)) to:\n  D:\\code\\bioelectrodynamics\\code\\NN_chemical_CAR_T\\Step_3_Programming\\enzyme_based\\brenda_ec_index.mat\nEC 1\/20: 6.3.5.1\nEC 2\/20: 1.1.1.1\nEC 3\/20: 3.1.2.6\nEC 4\/20: 1.8.4.4\nEC 5\/20: 1.1.1.78\nEC 6\/20: 1.2.1.49\nEC 7\/20: 2.6.1.15\nEC 8\/20: 3.5.1.3\nEC 9\/20: 2.6.1.7\nEC 10\/20: 2.6.1.1\nEC 11\/20: 4.1.3.7\nEC 12\/20: 2.3.3.8\nEC 13\/20: 1.14.11.2\nEC 14\/20: 1.3.5.1\nEC 15\/20: 4.2.1.2\nEC 16\/20: 1.1.1.37\nEC 17\/20: 4.1.1.21\nEC 18\/20: 2.3.2.2\nEC 19\/20: 2.4.1.11\nEC 20\/20: 2.4.1.18\n\nCoverage: K_M 23\/27, k_cat 12\/27, both 12\/27.\n  r_adapter_1    EC              missing K_M and k_cat  (no EC code in the motif file)\n  r_adapter_2    EC              missing K_M and k_cat  (no EC code in the motif file)\n  r_adapter_4    EC 1.1.1.1      missing k_cat          (see the Notes column)\n  r_D1_Ff1       EC 1.1.1.1      missing k_cat          (see the Notes column)\n  r_D1_Ff2       EC 3.1.2.6      missing k_cat          (see the Notes column)\n  r_D1_Fd        EC 1.8.4.4      missing k_cat          (see the Notes column)\n  r_D1_Fr        EC              missing K_M and k_cat  (no EC code in the motif file)\n  r_D1_Sf1       EC 1.1.1.78     missing k_cat          (see the Notes column)\n  r_D1_Sf2       EC 1.2.1.49     missing k_cat          (see the Notes column)\n  r_D2_F         EC 4.1.3.7      missing K_M and k_cat  (primary substrate unresolved)\n  r_D2_Fr        EC 2.3.3.8      missing k_cat          (see the Notes column)\n  r_D2_lk        EC 4.1.1.21     missing k_cat          (see the Notes column)\n  r_I1           EC 2.3.2.2      missing k_cat          (see the Notes column)\n  r_Nw           EC 2.4.1.11     missing k_cat          (see the Notes column)\n  r_Nr           EC 2.4.1.18     missing k_cat          (see the Notes column)\nReference snapshot applied to 0 of 27 rows (matched on reaction ID, equation and direction).\nWorkbook written: kinetics_BRENDA.xlsx\n","truncated":false}}
%---
%[output:812f1865]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Unrecognized function or variable 'row'."}}
%---
%[output:49e91d33]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAfMAAAEsCAYAAADXZ1qxAAAAAXNSR0IArs4c6QAAIABJREFUeF7snQW4HsXVxydQrAGCQ4DgQRJcg7QUt+JaKO4EhxSHIMEhQHB3C1a0QIqGj+LuHqQBgpciQb7nN+3czt2snN3Zffe9e888T55A3tmRc2bmP0enx6+\/\/vqr0aIUUAooBZQCSgGlQJelQA8F8y7LOx24UkApoBRQCigFLAUUzHUhKAWUAkoBpYBSoItTQMG8izNQh68UUAooBZQCSgEFc10DSgGlgFJAKaAU6OIUUDDv4gzU4SsFlAJKAaWAUkDBXNeAUkApoBRQCigFujgFOoH5Tz\/9ZP7xj3+YK6+80jz11FPms88+MxNMMIFZeOGFzeqrr27WXXddM9VUUwVN+fPPPzfbb7+9ee6558yMM85oLrnkEtO3b9+gNvN8\/MYbb5htt93WfPTRR2ahhRYyF110UdCcfvnlFzNixAgz7bTTmkUWWSTPUMR1fZqFjvmHH34w99xzj7nhhhvMs88+a7755psOHm+00UZmzTXXND179hSPreyK77zzjnn00UfN5ptvXnbTLWnvuuuuMwcddJCorzrWv2hgXbwSa\/rmm282d9xxh13jY8eONVNPPbVZbLHFzCabbGKWXXZZM9FEE3XxWVY3\/Fbvwbj+\/H208847mwMOOKC6CdfYMph71llnmdNPP90MHjzYbLXVVoVH0wHmo0ePNocffrgFpqTChjjuuOPMSiutZHr06FGo0yaBOTQ75phjzJ133mmuvvpqM2DAgEI0yfqoLDD\/9NNPzSGHHJLK4yWXXNKcccYZZrrppssaVqm\/c8m47LLLzCmnnGK22267Lrt5FcxLXRa5GiNlxmOPPWb2339\/e1lPKuzTE044wfTp0ydX+02v3Oo9mNZfdwHzF1980Qq3CM2hgq0F808++cTsueee5vHHH89cr5NNNpm59NJLC0uhTQJzDoTzzjvP0qzdwZyD7swzzzRDhw7N5PGf\/vQne7FrpfSCRshJ4135Jq5gnrm8KqswcuRIM3DgQKttyirLLbeclYamnHLKrKrd5vdW78G0\/roDmHOZOeqoo8w111xj0IoeffTRQWduj7Fjx\/7KokbUp6D6Q3r73e9+ZyaddFLzr3\/9yzz88MNmyJAhHbfdkI4VzPOfDWVI5l988YXZaaedrPnkN7\/5jeXxhhtuaHn83XffmVtvvdUce+yx9iDEZMCFbb755ss\/2IJftPogKTjMzM\/8Qwgt1qabbpr5jVYIpwBap91339088cQTtrEVVljB7Lfffmbuuec2448\/vvn444\/toXnOOecYVJsUzr211147vPOGtNDqPdjq\/tqNTa+88orZZpttDGu3jLXY45133vkVMR+7BUCODXmeeebpNG\/AZNiwYdbWirSGTRVAcAW7MXW+\/fZb8\/XXXxv+HwkeW\/sOO+xgbVTjjTeerZ4F5gALti5Uri+99JLtZ\/7557eTXm211TrdXNIk46SFkmYzZ\/w33nijueWWWwzqDzZ93Dz8NqILxD\/AuXnddttt5uKLLzavvvqqbYtL0i677GL69+8\/jqkCuj344IPm1FNPtXPHrIHKmQNnjz32sH4GRW3mPt0XWGABc+GFF1rQdoWxHnbYYdaWTuF3DkRfmmcNrLXWWp2m7C9Ixnn88cebSSaZxK6Fa6+91tLytdde66AlKs4tt9zSLLPMMnZNwO8DDzzQ0ila\/PYcbRgDvMHMs9RSS5kdd9yx0\/qiDR9Q0ZgwHkwHSG6UVVdd1QwaNMiu97vuusucf\/75neiN3aqo30BRMPf5w7zhxV\/\/+ler8WFvurXgjy1tHTpaMkcOCswX7Ik55pjDqvOiKubLL7\/c2uwobg1XQce8e0J6ALN+9tprL1udNcocuKj6hf3Mev7www\/tRXbRRRc1E044Yac6HKzMG5s7dE\/bsz790SaxV6Ej6x5\/o9lnn93suuuudv9GtVxoyp588knLi7\/\/\/e\/Wrj\/TTDOZ9dZbz9pN\/b3JAPOs\/7zjkuxBLvvOF4T9ghaXM3riiSe2ZzxnGmMEI4YPH27NHczJ+VxtscUWHed33v7iNHXO94f9gYBCScIJfnNYwX7gPOZyd\/bZZ9tvOUvcmTDLLLN0rAd4xDmM9pW9A08p8847r\/Ut4qIe5ZN0vfqa0jTfGZ+XaJNYv5NPPrnthjOWsxO85Hzocfvtt\/8KUFBgCIecbw9\/\/fXX7YJkYWcVl+bd\/x4wRqpHmuff08A8y6bLOPnjLhJlgjkaiCOOOMJu4rjiz+PNN9\/scKKL1nUHIZIwqmoIHS0scOqtv\/76HbR2jhBoSJzk4L5DQh4zZoy9wRUF8+gGWnnllQ2bBGCPHmj+eAFOHAZZyDgPoRby6\/sg4G6XeWj5\/fffZ4I548E8wAUjWuALqlX+uHXhgxD+HVyQojRlA80555xW6xQtIWaGMsAcfkNjLm\/RssEGG1h1HBcUKZgDGMwT\/w5K9FIGv7AzcxADQFzoZ5tttk6XojLomHdPZJ037vcff\/zR7rXrr7\/ergHWZBH\/FQAIaT7O3k67nI3sBbfOfPrjXPfVV18ZzoZo4ZLhr0\/WIjw56aSTxlmXfBv1W2Hv5ln\/eccF6GZdqH0w50L49ttvd0yTfb\/GGmvYdcWfpOLO77z9RcE8Cyc428AcH2gdVgB6\/fr1s0620TOBtc9FhXOBkmW2CfEv4jxHU4qDZhSkffoxBucUhxDEJd+tv3\/+85\/W3o6gyDne4+STT\/7Vqdijm5zJQgRuYBtvvLH585\/\/bKaffvpOvLrvvvs6bmwsCDqkswceeMDelLk9+INNAnP6cup+vj\/00EPtzeff\/\/63PbyQ8Pj3Cy64wCy\/\/PJ2DGWC+d\/+9jfrN0ABsLh8IDkmzSOtf99D0Z8LWgvaBuCjqmxuiGgfoBeL6sQTT7R+Ce+++651BnO3z6Jgzni5XbIxfJsi40MTw2GN5oPIAl\/r4h\/0Uaku6Tfmx8ZFuoVH3HpRdXIAABqAlA8abmxJNvObbrrJHjbQFSnASWCsFwA+ui58QHUmBS4irD0Oa+cbwm+O1xyAqGmdFMxhy00\/b8ljM\/e1OP6+oE\/4jKTMJoX3++yzjwUZDqO0sUUv31wYjzzySPP+++93qPSilzJfu+L\/ViYdi+wJKe19E1KS5iGrLaR1DldogTTu1i3+RAADa5r1guaJCxUleplCk4VAMPPMM5u7777b\/OUvf7H7GUmOC1Lv3r07QAIBid\/QInARwXaPVoz\/hlYc0ux7+sy7\/ouOS2rDZl9zPnFefPnll1brwCUG0OFs4eKy2267WakdyRYhkbUb5Y20Px\/MoQ3mQMyAFHcecFZzmYNX1PEvvdGzGv4efPDBVgvCmQxAcpGlgDtoWPh3zgPAlPOLvTjDDDNYTSLmGncxZl9SL69DuC8kRUHaX6toeRgrhXlvttlmHT\/7baD96XH88cf\/muTE5TYJUoxTn2ZtCvc7tycI\/cILL3SSJpPAHNBy6n5uVqianZosSTosE8yT5pU0jzQw5+DkBg94sRBOPvnk2Lm4hfDzzz\/bw4MNT4naT5AsUSezSEPAHM0JNkUOHFTfcQX1PwuVflzxpW\/\/wueDgH\/4OL5wELHp\/\/jHP3a6IMT1m7Sx\/U0FsLFW2VQU\/2bq+3H4IBRdS6ib3eXV\/8aX7mi7qENjWWB+7rnn2rVDge9caK+44gr7\/\/CDC3K0sF+56Dhzgi9ppF3KHH\/TLkWhdMy7J\/IcjmX4lPhrnP0BMLkx+GcTEj\/rB\/D1QZNLFhdLt2+g97777msjRwBxLmDY731zVvRi5n\/jgK9Xr14doCJd\/0XGxVqSgmv0TEvDBF8VHFUnS\/vzwdw\/c1iT7GfAmeLTlgsHgI\/GJHpWI5RyaXJCCwIpWEXhb8wJvvYV4ZEzETNIWcUJPLTHOuFCEC1Z+95vgwtUKpi7TcIClsT5cWtBMkCS4GbqnFF8AEoCc1+dIPVmrgrMJfNIA\/O0RerP3x2QjqFs\/Dj7iQ9aIWDuFguLhA1\/\/\/3328PG+Qe436MRC\/5h5ktuEhCgTQ4zNtA666xj5xd3UEt8HHwbOu36pgNf+vEB1W1O16f\/W3QTlRGdUAaY+4e\/44k\/tjjHOl\/y5RsOqtNOO83a9lzx7cruUuaDfBQsyqRj3j0RtXenHaKhYB61PUY1H\/5FD5Bg3aM180Ez6ocSNWu5y2GacBA3R78P6fovMq48YJ51PkMvzgw0cGg8MfGwPssAc19SdVK0T7ckkPT3T1TC9demozFrwqnB3X76wx\/+YDXPiy++uDVzhRSJOc4XZPBx4XLCBd0VXzBhP3cCc35ELeeKW5BTTDFFJ129PwmkvUceecTaShx4RycpAXN\/clFiJxGtTDDPO480MJce6I4utOUS6cQ5p4UeVlmLjs3HJQwnMZdnwAdt\/8brJAYkE2dn9aUV+uLAQkJkk0TLEkssYbUMOAL6TkFJh73\/72nz8A+KtLCWtE1UNpjn8WbPcgzNAvOofQ\/VJpcV32QSdyl76623OtTvUZVhmXTMuyfyJKfyaRd36GWtf\/\/7qEo87jLlLkJpzrRJYO5\/g4odjWeas2WR9V9kXHnAPOl8fvnll+0FEtNk1B5N+2WAub8P0ISsuOKKndjr08s\/w9L2TxyYA9a+ecPvBGEHkyBOfZhUihQfiJPOCV+bhYYBe74LpeRMRv2Oj5e7YPa46aabfsWO6ADF2WncAEmIgh0Czz9f9ep+9w8RHLvIFIczBCol2o16YCcdWv5mj7txxRGsTDDPO48ywNwBN\/ZkB+Zxh0kZYC5ZPGkaAP\/Gy2EGqLuwiji7kfNmd97YUf5FnVRCwdyXZssEobwbVXLjjmszBMxJXoSPgvOrSHLM8SVMdylDO4M9MKqWZIxl0lEK5nGX2Swe+Adb1FSQ9S2\/1wXmUdNF3FilYO6v\/6rBPA58fF8NeMBlncsKGgxMppwfXQ3MnTc72OciDqKgjvodM2IesxBtSM5jhGMuDFyMon4uvobHYUaPl19++Vd3KEdD0wg1wIsZWy4pXrmR4PRASAc3BA4HfsMzEGcHJAFsTUhceW3mvpo9qhqVSOZRrYIPPr5aKG6hcwtz3pwcarj\/Oye7IjbzvCYDX80Xp2ItQ83u0yNqa3L09SW3qDofwICO+ECwsGaddVbrjeurHeP4xIYgDASgYUPwxzng+ekLk8CcQwL\/A2gQVTMmrYsyQUgCBn6dVoN51CGIPYpNN86mzjh9\/wskLOyFaGPiPGrLpGPePZGX7r7NOyk0jbWI3wFrjYgFgAZ1fhlq9uh+kajZJReXIuu\/1WAOXTkL8POgAPacEQBc2TbzMtTs0ctIkmTur0GwDmc+sI5wVv6bEjVNSdet5JxA6naCdtQk+Mwzz1isdU6UaHhs0hjfO5BDGu8+Njf2NIAem6qk4BV91VVX2VznvtOaRM2e5GRCv74bvy+5+rcbX\/UTzXaWBeZOKxGXLz5pHmmSeZqzj+9s4UIBacvfDFEHOL6hbogD3AcffGBtQIQxUDjM8DjlAge9oD8aGCQ1StTJxUU2cHkDwLm8cdHL4wxDu3iNMheKz5ciDnBRb07npFkmCEnWfZ1gHlUFRp23ouP3L2U4YBEayIEQlxe6TDrm3RN5JR3fG505c2Hlgk6IHW1xoWQ+SFmERlH8OXOBR4KkFHGAk4J5mpNWFBRRI2OfdV7VUeBIWv+tBvOkiwu09AWRMiRzH8SKOsAVAXN\/H2HLBlw5l4u+r+DjQJwDHGuBcbpwXL9OEr6lpnPlIzaCix9nQkjd3HyxFXDAc6vn9khBagcQOBwIh3HA4N9AJaFptOXCG\/hvPxSAf8ceiyrHlzZRLeI1DnEBDMwFTgLMAnNfMqc\/LjPYRDgAkuZBPf8ywQGAkxeHIzRy3ul+aBSXIxcGQWgawOjCn1ikXJwYsx+ahtc5thEXc1zUAQ4eIr0wn6wS5zzFNwAut0HfHhYFgehhBV1wfCNumrkRL+vCSvxLi69ScpuUcBN4g1bIJTRBa4NKmd+gP20xXqQCZz8rE4SyaBX9XapOdt85x6giavakMLQ057EoYDCOpIx\/ZdLRvwxK90Re2mfFBfvt4efB+nOxyFxmuezixZwWmuY7FRYFTcbpQtP8c4s94OLc8S3hgsH4fK2DZP0XHVfaHiQRjEsaEwXD6JriHGMehBVDL85vSlTrKO3PP7\/9FKi0mRSaFs0Vkddmjv3fhQqjpSWaBG92FxHk+BT1F5KuWV\/jEhea5ms0aBNMYL2Sg+Lee+\/thG9OK93x0Ao3ZwAwzmnJDZCbPC79vo0gaQNhP8dxDjW1H4IRkjQGbQGA7R4BSXO0IsMYQIj0mAXmaBLyzgOa+CpLRyPXF\/GpOIi5MCGfydFEOvzGIiHOlDSrUecRNjaSR+hLb2nJJ9z44pKwuN8IfeIy5dZINFbc1SMkz\/lLJC1uF\/\/sgMeXGN037uLCRQDbFHGkcSWaTKhMEJJuTlevlWAu7St6+Po5DRh3UnrmsumYd0\/kpT17CDMOoOOydcW1wX6CJvgN+CUkaYxUMqe\/rKQxnJdcep2pBCEgz\/ovCuZpexAASQJz5pSUXIzzhAsJEjrFd1qT9ldm0hjGIJHM3cMnSYl9aCfLpJW2fn1P9TgTl6\/R8DWhtIkzHOG5CLO+mbPTE6jcevBMZxNzYPvPY+LYhuQNgDqJF1BFQn\/66aetnY5DAknKpSxl87r4Pae6AhDSnkB16VyRWJx6H+mVm1bc85yMh3zLpGHlW5fSb6655rI3bQBQAubOm106DzQWbEoyI\/GNSzbCwYhDGADknMAYG+ptFgiMY1zcyKOqxGg6V26CaAjKSOfqFhZ9cCNmTGhO3KEHMBN6ga0L58UkNaevjkzL0e9S4\/rPUDJ\/UrCSgAj1fDTFJfZ4Ng9rj\/659WI+gQ6sTcIdkcT9dYHKHonc99ouG4TygIoUYF2bIZK5tK\/o4eUfJIwjKS90FXTMuyfy0N7VRWAgFS6heC7sUvoEKt9yqSZJVVYK5qKg6S7vpHMlIRehW5yprHP2Ons+6iWdZ\/2HjCtpDyKUpIE5c0Ka5XzgQkVx6ZYR6kjEwlmDtIum0e1XSX9J6VyhG1oLLmGc35J0rlIwB8ecAxx9YCtHoKLAJ86vrbfe2vjpX\/OuVXeWxkVg+CYUBB8S3CBIkwcEPuCEjHDjJ+LpBOZJg2EDEntJA87e5OpyQKN6RTpyOWPzTkrrdw0KRG01cbnau8ZMuvcofQ1L0Yxp3ZuCOnulQDgF\/AQ40TC7aEIY5wiX1msmmPt5tlErkPbTPdKAat55J0fVpuFT1RbaiQLuloq\/AlqIol6c7TSn7jgWtElJaUO7Iz10zkqBuijg2\/+jdnNfAyoVmjLB3GWNco4G0QQHSO0uR3aRZ9xGjRplPUxx90dNjj0AL2v\/pTVHbPpC5cGrN6hsUDng5YeDgnuVzdV1rzOhggd80l4wqouZXaVf33mEMcfZ\/LvKXLrrOKO536EDe4L9gVlFi1JAKdB6CjjH52mmmabjgaO4hDCEUWaVVDB3IQfYpV0+4rgGndoOb0VsnFFbaNIgnNMEdlwuC3x\/++232z\/+S2t87zQEDz30kLVV8PJNUl2kD24zeDgTo4x9GpsH9im8SP2X17IIpL8baz\/k0QjoygUKmz+JgXw7tdKpvSnAAYHfCnY2TGNo2Ah5Yh9pUQooBeqhAGcqIZFglYsMkiQxihttKpjnyc2O9IbjEuFWklSMLoQJpzlfOnCqB9cWcaIUnHHw6PQTujhC4Bjlnm2krgv9IJOcS16AmpiDjDb8BPz1sFB7VQooBZQCSgGlQHkUqA3M3Zu0eG363o0OuPHYc56+7jEIJHj\/NTXqOq8\/pG0XA83FAq98\/9lB6jp3f1T40bS15ZFUW1IKKAWUAkoBpUBrKZAK5i6XM48xoLZ2T09Gh+hyQ\/OoO+9DE5ZVtDiPaWzj7vUilz2KECRA3g+bctoDss\/RN0lbUB\/ydGD0EQNnNsAsgISv3vdFuaTfKQWUAkoBpUA7USDTAc6lEiXGnIQmUUAHyLFvo+qOPrmYd6JI4LjkY3fHOx7gxv7upG\/c8\/3H2WnfZcr56quvLEDz\/8SxE88dBX6XIo85cVFwXvl5x6n1lQJKAaWAUkAp0E4UyARz50yGZE7BA9ZJtCSfwFOcEs3ClWeS0UfYsXNj73bZwVze7rjXepy0jVc8anW83HF6I\/lC3BvsqODx0AfM+\/btm2eYWlcpoBRQCigFlAJtSYFMMGfULvsQjmqkSPUL6m08xJHcpV7sUUqgzn\/88cfNzz\/\/bL3OAVqyB7nUrXWCOQ+U8EeLUkApoBRQCjSXAvhvFX2fvB2oIgJzN1DU1KQddLnDCU0ikUzeF47SJk4f5AHG7k0IGtI1aQ9JcdhqyRwQHzRokE0XqEUpoBRQCigFmksBBEjSSXdVQM8F5q1io8sdzcUB1Tn51VGdp9nM0R7g6c43WTZzQteinu5xc3MaARhMPt5WFS4PJOCpqt+Q9slPzNjK1Fawecgsl4fGIXOQ8rHqPqpun3k2oQ+dg2zFVk2nqtuvc726ubkIKhnF26vWOA+tPP\/883aECy64oP2b\/wcoJQU1O98VVbe7PqJ2cJzaAPNVVlllHDt4kjf7lFNOab3ZSZgfbVfqze7AvNUMrrrfkPbdtz3XGmjG6zVt4rL4adRL5ruRw420Xl4ah8xBspapU3UfVbffijm0oo8m0EnnINt1ddGpFf3KKFC8Vicwd8BIc0iuFPfCmaSLPG9t8z43Kmy803lRxy9jxoyxWcZIcYe0TeE50VbHmdfF4Kr7DWnffTvZ5oPNBLP0T1wWY0e9ZL65erCR1lMwHyDZYrnrhPBa2lnVfVTdvl5IZJxuAh+SeN2KucmoXLxWJzAnNMy9GY1HOYX\/J+xLUojt5jvnhZ72jUsByytsfgy7y9RGWJlvI0\/KAEcYG++KSzLAEbuOU500A1xdDH733Xdt6ll8BlwGPAn9pXVC2peC+S9ffWp+eOF+M9ECK6RK8A7084J5yBxaQSdJH02YA\/Oseh5Vt69zkKzW6vncCj4k9VHXWS+jvKxWrTZzbNcDBw606V9dbvY777zTvhG+wQYbdApPc7nZR4wYYXbccUdRbna87wHEpZdeulBu9roYTOIbstWRq37iiSeWcTJHrZD2pWAuHU5RMA+Zg3RsVfdRdfvMswl96BxkK7ZqOlXdfp3rta6zXsZZWa1UMHeS+owzzmhWX3311BaRJEmhip1akpvdNcZjKyR7QbrGeY2HPADruFA3HON4Me2qq64SvZo2fPhwG+ZW9NW0uhhc9aYJaV\/BXLaxJLVC+CBpv87DUTo+Sb0m0EnnIOF0fZfPus56GVVktUrJze4yqz3xxBPih1Zkw6u3Vl0Mrnrjh7SvYF7emgzhg3QUTehD5yDjdtV0qrr9Oi+fdZ31Ms7Kao3jzX7KKacYQJlCMpc333zTZnxLi71z9XhXPPoQimwY7VmrLgZXvWlC2lcwL2+thvBBOoom9KFzkHG7ajpV3b6CuYzPSbXGkczJg44HO6+a5SlTTz21GTp0qM2J3pSiYD4uJxXMy1vdTT4cy6NSfapXnUNnCjR5vdZ11pe5xsYBcz\/LG17se+65p1l88cVt7vXEG0GPHjYTHBnhmlTqYnDVmyakfQXz8lZ4CB+ko2hCHzoHGberplPV7atkLuOzWDL3K+ZxgAsbRnt+rWCuknm7RhVId0wTDmCdg4zbVdOp6vYVzGV8LgTmYU13\/a8VzBXMFcyz93HVh3zV7bcCRHQO2euoFXxI6qOus15GFVktUZw56VVfffVV+1Z4tJD8BSc5nODI133kkUfmCk2TDbOeWnUxuOqNH9K+qtnLW4shfJCOogl96Bxk3K6aTlW3r2Au43NhyfyTTz6xdnOeKM0qedK5ZrXVDr8rmKtkrpJ59k6s+pCvuv1WgIjOIXsdtYIP3VYyxxmOhC6kW11nnXUMoWdnn322mWuuucxGG21kXnnlFXP\/\/ffbUDbSpK6xxhqdHjaRsa99aymYK5grmGfvz6qBqur2WwEiOofsddQKPnRbMMcBbt9997VcIP6cF8iOPvpom2rUxZPz+AmZ1s4\/\/3z7Z84555RxrQvUUjBXMFcwz96oVQNV1e23AkR0DtnrqBV86LZg7l5RGzBgQMfTozx4csEFF9g0qX369LEcIg0rb40vsMACZvfddzc9evSQca7NaymYK5grmGdv0qqBqur2WwEiOofsddQKPiiYe2AOwJE7nRfIFllkkQ4OnXDCCeatt97SDHCyNZtaq+qNH9K+OsCVwOD\/NhHCB+komtCHzkHG7arpVHX7CuYyPifVSvVmx3v9wAMPND179jRHHXWUmXDCCQ0Po2y77bbm4IMPNmuttVYnMOeg5ynSPA+thA2\/2q9VMlfJXCXz7D1W9SFfdfutABGdQ\/Y6agUfuq1kzsSRwLGX8xraKqusYr788kuz0047mZlmmsm+N87b5aNHj7YZ4rCp4zBHLvcmFAVzBXMF8+ydXDVQVd1+K0BE55C9jlrBh24N5n5o2jLLLGM926+99lpz0kknmXnmmccstthi1qOdQ2\/QoEFml112UZu5bN0m1qp644e0r2r2QOZ6n4fwQTqKJvShc5Bxu2o6Vd2+grmMz4XU7O4jJG9AHGYec8wx5ueff7Z\/X3\/99R3trr\/++ubQQw81U045ZdiI2uhrlcxVMlfJPHtDVn3IV91+K0BE55C9jlrBh24tmSexgJC0Dz\/80Lz33ntm+umntyFp4403noxjXaSWgrmCuYJ59matGqiqbr8VIKJzyF5HreBDtwVzksYgfc8444xm2WWXbRxYZy0vBXMFcwXzrF1S\/ROlCoTZPGgFEDaBD90WzL\/++msbNz777LObww47rHFPnGZtEQVzBXMF86xdomCeTaHqaaRgLuHCf+rEXUrqOuvlo86umRqaxiMqeK7znvkBBxyQ3VrDatTF4KpvwCHtqwNceYs8hA9HsdM5AAAgAElEQVTSUTShD52DjNtV06nq9ltxIem2YI6a\/dxzzzU33HCDDTmbb775ZKuqIbUUzFUyV8k8ezNXfchX3X4rQETnkL2OWsGHbgvmP\/zwg31EhZzrANv8889vllpqKZs8Jq706tXLbLLJJjb2vAlFwVzBXME8eydXDVRVt98KENE5ZK+jVvCh24K5y83+3HPPiTihT6CKyJRZqeqNH9K+qtkz2SeuEMIHaSdN6EPnION21XSqun0Fcxmfk2ql2syRzJ9\/\/nnD35Iy0UQTmQUXXNDwdxOKSuYqmatknr2Tqz7kq26\/FSCic8heR63gQ7eVzGXk\/1+tn376yYavNSXeXMFcwVzBPPsUqBqoqm6\/FSCic8heR63gQ7cFc7zZDz\/8cNO\/f3+bpjWtkL995MiR+tCKbM2m1qp644e0r2r2Ehj83yZC+CAdRRP60DnIuF01napuX8FcxudCava498zjGkIi51W1F154QcE8jB\/266o3TUj7CuYlMFjBPBcRQ9artKOq+6i6\/XY\/N9qFD91GMgeUjz32WHPppZdKad+p3tprr21fV+P1NGnhSVXC3h588EHzzTff2NfYNt98c7PVVlvZp1f9wpOsvOJ28cUXm88++8xqDPbdd1+z\/PLLj6Pax85\/2223mXPOOce88847NvHNrrvuahij1KavavZxuahgLl3Z2fWacMAriGTzuRU0akUfTV6vdZ31stUjqzWOA9y7775rE8W8+eabshaMsZnhBgwYYAYPHmzmmGMO8Xeo5QcOHGjfP99hhx3MzDPPbB599FFzySWX2PZOPvlkM+2009r2\/vWvf5kjjjjCPPTQQ2brrbc2\/fr1M7fffrv9M2TIELPRRht1vNbGpYSHYYiR5+31JZdc0rZ72WWXWUDnuVbGnFXqYnDVmyakfQXzrFUj\/z2ED9JemtCHzkHG7arpVHX7rbiQdBvJPLpkpGp22VLrXAspnFSx\/\/73vy3wzjDDDB0V7rvvPmuj5xU2JHTKddddZ19qQ4pHEqcA2qeeeqq5++67rXp\/ttlms\/\/OJQHQ5nvi3nv06GFcnnnaQPPA061ZRcF8XAopmGetGvnvTT4c5VTIrtkEOukcsvmsYC6jUVItUWja5JNPbt8uL7Og+t57770Nb6RHU8W6S8Qss8xi1fY8ubr\/\/vsbXmoDvP2kNC+++KKVvpG2AX4A\/oQTTjCPPPKIBfjevXt3DBvP5O23394+GkOfWdK5grmCuXqzZ+\/6qoGq6vZbASI6h+x11Ao+dFvJXEb+8mtFwXzMmDEWsFdccUVz0EEHdajT6dnV5bKBEx6bBomfbHRcBHy7+3fffWcOPPBAg5c+Ej6XlLSiYK5grmCevb+rBqqq228FiOgcstdRK\/jQ7cEcwLzrrrvMW2+9ZT3WkX7jCjZvAHXKKaeUcS6h1lNPPWW22WYbs+WWW5pBgwaZl156yYL5fvvtZzbbbLNOX+EUB0B\/9dVXFqD5f6Tv5ZZbbhzgR9V+3HHHGdT42OX79OmjYD7xxLl4pWr2XORKrdyEA74VB3AT6KRzkO2buuhUl+Amo4qsVqqanSbwNsf+jFo8q5SRzhVHN6RvHN2wbS+yyCI2Lzwe7gDxpptu2mkYTtoeNWqUVavj5Q7w47Ue99IbKni83AHzvn37KpgrmCeugaoPlqrbbwXQtqKPJtBJ55CFHv\/5vS46NR7Mnf0ZL\/CNN97Y\/PnPfzbTTz99IldwNEO9XTQDHOFkSNdnnXWWlchxgqPNusF8r732sg\/MUHDU8531ZEs0Xy0W9OjRoy2t84T5SXsJaf\/JJ5+0vgmTbT7YTDBLf2mXifXGjnrJfHP1YBtuSASDtITMoV36aMIc3AHcrutVeS2lQHa9pq1XZvzll1\/aiX\/44YcWc66++upc51A21VpXQ\/Se+Ywzzpg7fjzvFJCwhw4dai688ELrzOaHj9UN5v5cCItzHvZ55yitz6Xm008\/tWF50ph4advUC2mfR3eI7S8bzHFsRLMjLSFzaJc+mjCH0PUk4UUT6KRzkHA67GyS9fC\/Pu68804bJeWXxoJ5laFpPgG\/\/vprc\/TRR5u\/\/vWvZrfddrNqfR\/EnMd6ms2czQIgEPKWZTN3aWd9T\/e4heAuESeddJJNZkNphWTOxebjjz+2fU2cUw0uWdAh7UOT7bbbrnQwJxnQ4osvLhm+rRMyB2knVfdRdftKJymnq19PTeB1E+bg7wlfMn\/sscfM6aef3lzJ3Nmjp5hiCnPYYYdlhnLJt87\/aiKBHnLIIeaBBx6wceHYxqMhY++\/\/761g6+yyiqJYWxRb3ac8KLZ6NSb\/T90D7FLqQNckVUe\/00IH6SjaEIfOgcZt6umU9Xth55NMirFn3+Nt5lDHFQRZFg7++yzc6lBJYT95JNPbPw4tyKSuWy44Yax9nac4jTOXEJRWZ2QTalgLqOxpFYIHyTt13k4SscnqdcEOukcJJwOEzRkPXRTMCdJC97heIlfc801FmwJ+Zpwwglj6ZbnPXOXB\/7KK68cJx1rXONJGeDIJU9ed0kGOFS5pIjVDHDfm6Lx0wrm0iMju14TDvhWXBiaQCedQ\/Z+aMVaSuqj8ZK5s5nj9CQpeULTaBNHMmzC6667bqzXNnHg66yzjr08uNzsI0aMMDvuuKMoNzuPrOCwtvTSS2tudo+BIQeLgrlkJ8jqhPBB1kN9ko50fJJ6TaCTzkHC6frWa+PBHAC9\/vrrbUIWSSEsjVzofrrVpO+QtIknTyvRV9hwlCOE6aqrrhK9mjZ8+HAbT66vpnWmcsjBomAu2QmyOiF8kPVQ3+EoHZ+kXhPopHOQcLq+9dp4MJeRv7m16mJw1Rs\/pH0F8\/LWewgfpKNoQh86Bxm3q6ZT1e2rml3G56RamRng4j7EKxzGhiSICRt2a75WMB+Xzgrm5a29Jh+O5VGpPmlN51CeRk9Ky7r2RF1nvZQuknpiMCddKh7t5GgnltvZx8lQRrKXtdZaq3Rvd8kEqqxTF4OrXtAh7SuYl7fiQvggHUUT+tA5yLhdNZ2qbl8lcxmfC0vmPE5y8803W\/v22LFjO9pxYM4\/kKQFe\/b5559v5pxzzrARtdHXCuYqmRf1+pcs4yYfjpL5S+s0gU46Bxm366JTXWe9jCqyWpmSOZPceeedzbzzzmv+8pe\/WOn7lFNOsfnSCQcjOQuPouy5555mhx12sM+Pkk+9CaUuBle9oEPaV8m8vJUdwgfpKJrQh85Bxu2q6VR1+yqZy\/hcSDL\/8ccfzeGHH26fPUXqdilNeXnMgflUU01lc32TIY6k9aRUlXizhw27NV8rmKtkrpJ59l6r+pCvuv1WgIjOIXsdtYIPSX3UddbLqCKrlSqZJ+Vmj4I5XcX9m2wI7VurLgZXvfFD2lfJvLz1GsIH6Sia0IfOQcbtqulUdfsK5jI+F5LMFcz\/8456q1\/SqXrThLSvYB624fyvQ\/ggHUUT+tA5yLhdNZ2qbl\/BXMbnQmCO+vzggw82eLIPGzas4x3vqBSO8xvPYhKqRnrVKp7tDJtmsa9VMlc1u6rZs\/dO1Yd81e23AkR0DtnrqBV86LZqdibOQyt77723+cMf\/mCOOOIIazf3wbxnz57mzDPPNGeddZZ9Qo6sbU0pCuYK5grm2bu5aqCquv1WgIjOIXsdtYIP3RrMSenKi2akdeVpUkLPxowZY183W2mllewjJxx466+\/vjnyyCMb4\/wG0xXMFcwVzLMP4aqBKrT9Dz74wPAnrfDwE+cazz3zXkRamXnmmQ1\/8pTQOUj6qrqPqttXMJdwOblOZmgan6JuJ5c6qnZeUfPLZJNNZh8+4b1xpPQmFQVzBXMF8+wdXfUhH9I+ID5o0CD7zHJZZamlljInnXRSLkAPmYN03FX3UXX7CuZSTsfXE4G5+5Tb68cff2zee+89w\/Oo3E5nnHHGxCdRw4ZW\/9cK5grmCubZ+7DqQz6kfbeHe6410IzXa9rsyWTU+GnUS+a7kcNzO8WGzEE66Kr7qLp9BXMppwPBnHzso0ePNrPNNltHUhji0EkYM99883XEoIcNp72+VjBXMFcwz96TVR\/yIe1r9EU2\/6Q1QvjQ7n3UddZL6SKplymZk8719ttvN0cddZRZZpllzPHHH9\/x9jj2dLzYsZtjV99www2tLb0ppS4GV71pQtrXw7G81R3CB+komtBHyBx0vUpXSna9ED5kt\/6fGnX1UddZL6WLpF4mmLt0rhNOOKHZbrvtzDbbbNMB5tjS77vvPpsd7u233zbnnXeeGTBggKTfLlGnLgZXvaBD2tfDsbylG8IH6Sia0EfIHHS9SldKdr0QPmS3rmAupVFSvVQwd+lcn3jiCXPOOeeYueeeO7Yd1O977LGH9XRHggf4m1AUzFXNrmr27J1c9SEf0r6CeTb\/pDVC+NDufdR11kvpIqmXCuZffPGF2WmnncwiiyxiX01Le0CF2PMnn3zSSuk8vtKEUheDq940Ie3r4Vjeyg7hg3QUTegjZA66XqUrJbteCB+yW1fJXEqjQpJ5UjrXuMY0N3soK\/73fdWbJqR9PRy7Dp8ZaQivpTOtuo+Q9nW9SrmYXS+ED9mtK5hLaVQIzPFgP\/DAAw1\/8+wpMeVx5ZtvvjH77beftaX7DnKhg6v7e5XMx+WAHo7lrcomH47lUSnsQqLrtTxONHm91nXWl8cdYzId4G677Taz11572bfK+TuaGObbb7+1aVwvvPBCc9xxx5lNN920zPHV2lZdDK5604S0r4djeUsyhA\/SUTShj5A56HqVrpTseiF8yG5dJXMpjQpJ5nxE+Bk52W+++WYrmc8777xm0UUXte09\/fTT5tVXXzVI5prONZQV\/\/u+6k0T0r4ejl2Hz4w0hNfSmVbdR0j7ul6lXMyuF8KH7NYVzKU0KgzmfIj0fckll5gLLrjAArdfNJ1rKAvG\/b7qTRPSvh6O5fE7hA\/SUTShj5A56HqVrpTseiF8yG5dwVxKoyAwdx+TzvXLL78077\/\/vs3XTiganus8wNLEomr2cbmqh2N5K73Jh2N5VArTLuh6LY8TTV6vdZ315XFHYDMvs7Ou1lZdDK5604S0r4djeas4hA\/SUTShj5A56HqVrpTseiF8yG5dJXMpjUqTzFGzk+I1rhCH3qtXr8akdFUwV8lck8ZkHzFVH\/Ih7SuYZ\/NPWiOED+3eR11nvZQuknqZ3uyo1m+99VYbmsbBllYWWmghc9FFF5mppppK0nfb16mLwVVvmpD29XAsb9mG8EE6iib0ETIHXa\/SlZJdL4QP2a2rZC6lUWHJ\/M477zR77723AdRnmmkmM8000yT2yZOopHMtmgHu3XfftRnnaCMuxzuOeJdffrm5+OKL7bvq\/fv3tw+9LL\/88uNoA7DpE1ZHGtp33nnHzD777GbXXXc1a6+9tplooolEdFMwV8lcJfPsrVL1IR\/SvoJ5Nv+kNUL40O591HXWS+kiqZcqmbukMU899ZQZOnSoWXzxxVNTuko6TKrjh8BdffXV44C5+50nV7feemvTr18\/+5obf4YMGWI22mijjrFx8Rg2bJg599xzzbbbbmuWXHJJ8+ijj5rLLrvMAjp55CVOe3UxuOpNE9K+Ho4hq7zztyF8kI6iCX2EzEHXq3SlZNcL4UN26yqZS2lUSDJ36VyXWGKJzNzsIQP5+uuvzamnnmqlbkocmF933XX2mdUzzzzTSuIUQJvv7r77bqve5611ysiRIy1oH3rooWaTTTaxII+d\/\/rrr7dtXHrppWaxxRbLHLKCuUrmKplnbpPKY9lDQETBPJt\/0hohfGj3Puo666V0kdQTPbSCRH7AAQdI2stVB4B9\/vnnzbHHHmuTzyy88MLm4YcfHgfMkcr3339\/88svv1jwnnTSSTv6efHFF630jbS91VZbWYAnT\/wjjzxiAb53794ddTmYt99+e7Psssva+WRJ53UxuOpNE9K+Ho65lnhq5RA+SEfRhD5C5qDrVbpSsuuF8CG7dZXMpTQqJJk7YHzvvfdSc7MXHYST\/D\/99FNz8sknG\/rhdbaoZE5cO4C94oorjqMhcG3MM8881tbOgtt9992tVz154v30s85swGtwSPiTTz556tAVzFUyV8k8e3dXfciHtK9gns0\/aY0QPrR7H3Wd9VK6SOplerO\/\/vrrVmWNTXqbbbaxj6mUVQDVa6+91rY97bTTGlTpcWDupG8ec9lss806dY9THI\/BfPXVVxag+X+k7+WWW24c4EcTQP74++67z2a069Onj4L5xBPnYqcejrnIpZJ5CeQKARFdryUw4L9NhPBBOoq6+mg8mKPexs5MDna82ieYYAIz11xzmQknnDCWN6He7Elg7ggd95CLk7ZHjRpl1ep4uSPF47UeZxpABY+XO2Det29fBXMF88Q1UPXBUnX7TKwJfYTMQcFcCqPZ9UL4kN26qtmlNCqkZncq7Oeee07UT2icebuCOa\/FLbXUUpYGM8wwg\/1TZWHTjB492kw\/\/fSlakLcmEPaf\/LJJ61vwmSbDzYTzNI\/mAxjR71kvrl6sA03jAtHTOogZA7SQVfdR9XtOzCvci21oo8QOul6la727HohfMhu\/X9g3qr1So+kJ6d8+OGHZtCgQbHO19Kx110vVc2Owxnq66SMb9HBh2aAa1cw9+dJWBxgVmUhRh4\/AkwP0pj4POMJaZ+LHbH9ZYM5jo1cBqUlZA7t0kcT5gAtq55HSPu6XqWrPbteCB+yW\/9PjVb2gbYZzPFLXCSVdOx118u0mbdygCE2cxYBgEC62SybOaFrUU\/3uHk6Fd1JJ51kE+ZQWiGZYzr4+OOPbV8T51SDS\/gV0j402W677UoHc8ISiZqQlpA5tEsfTZgDtKx6HiHt63qVrvbseiF8yG79PzVa2Qf9Ocn8scceM6effnpzJXOfAUjpb731lk2+gl2ago2ccDJeTxtvvPGk\/EqslwTmzpt9lVVWGccOnuTNThY6vNl9hz31Zv8P6UNsX2qDDF7mHQ2E8EE6iib0ETIHXa\/SlZJdL4QP2a2Hn00hfTTeAc4RB3vCkUceaUaMGBFLL2ydOJZleYdnETsJzDXOPIty+X4P2ZR6OOajdVrtED5IR9GEPkLmoOtVulKy64XwIbt1BXMpjZLqZarZP\/nkE7Pnnnuaxx9\/3JAJbv3117cSOeWDDz4wN998s3niiSdsytQzzjjDTDfddIXHlATmNJiUAY6EMw8++KAoAxyqXOLZNQPc9\/bRHBLq5FXj6+FYeHmP82GTD8fyqKSaJHfWct6mFfKCjBkzxkwxxRSZ+5oz3J3jUl41eb12C8kcACQF6tFHH2023njjcdTpqN+HDx9uDjvsMJs+NcQ5LA3MXW52tAM77rijKDc7j6zgsLb00ktrbnZvx4ZsSgVz6dGXXS+ED9mt1yvpSMcnqRdCpyasV0AcT2vsumUVonPwBcoD6CF8kI67rj4aD+YAKJ7L2MORaP00qj5znBqcOPRo1jUpE530HZc0xrVBDndCmK666irRq2lcMogn11fTOnMhZMM04XCUrskQOkn6qLp9xtCEPkLm0IT16ubQc62BZrxe00qWVroEP+ol893I4bmdvUL4IB10XX00Hsydcxk28azc7NjMIYi+Zy5dtsn1ql7QIe034XCUciiETpI+qm5fwdzYM2nzzTcvPfoibwhTCK+bMAfJfqhzvTYezEm3yvviiyyySOqraS5N6jPPPGPOP\/\/8wu+ZSxneqnp1MThk40toE9K+HiwSCsvqhPBB1oNK5k1Yr02YQ7uv17rOeildJPVSHeBwqMBWTlw29ue55547tk2Xv5186NjOs14jkwysHerUxeCqD\/mQ9vVgKW9lhvBBOoom9BEyhyas1ybMod3Xa11nvZQuknqZ3uxI2zywMscccxgeOsFr3eVm\/\/HHH62X+ymnnGLefvttc9555+VKySkZYJ116mJwyOEloVdI+3qwSCgsqxPCB1kPKpk3Yb02YQ7tvl7rOuuldJHUywRzVOg33HCDOeSQQ+xb4XEFSXzIkCH29TNSujal1MXgqg\/5kPb1YClvdYfwQTqKJvQRMocmrNcmzKHd12tdZ72ULpJ6mWBOIwA6qnQ8w++66y6bMpUy2WSTmTXWWMO+UoYKvklAzvzqYnDI4SVhekj7erBIKCyrE8IHWQ8qmTdhvTZhDu2+Xus666V0kdQTgbmkoSbWqYvBVR\/yIe3rwVLeSg\/hg3QUTegjZA5NWK9NmEO7r9e6znopXST1MsEc1XqcQxsZxIg7RzpvaqmLwSGHl4QXIe3rwSKhsKxOCB9kPahk3oT12oQ5tPt6reusl9JFUi8RzEeNGmUTwDjHN1+F7rzcr732WrPFFluYXXfd1T7X2bRSF4OrPuRD2teDpbxVHsIH6Sia0EfIHJqwXpswh3Zfr3Wd9VK6SOrFgjmhaPvss4\/Nsrbyyivbp0X97G\/ffvutdXi78cYbzdixY81cc81lzjzzzMTQNclA2rFOXQwOObwkdAxpXw8WCYVldUL4IOtBJfMmrNcmzKHd12tdZ72ULpJ644A5L6SRKObNN9+0Wd823XRT07Nnz9i2cITj0ZJhw4bZt6h5D7ZJEnpdDK76kA9pXw8WybaS1Qnhg6wHBfMmrNcmzKHd12tdZ72ULpJ644C5e1gFSZv3w7M81P3QtdCHViQDbmWduhhc9SEf0r4eLOWtwBA+SEfRhD5C5tCE9dqEObT7eq3rrJfSRVKvE5ijPj\/wwAOt6jztYZVow6R9HThwoH1SkydJJ5poIknfbV+nLgaHHF4Sooa0rweLhMKyOiF8kPWgknkT1msT5tDu67Wus15KF0m9TmCe52GVaOP60IqE3LI6VR\/yIe3rwSLjoaRWCB8k7VOnCX2EzKEJ67UJc2j39apg7nFIwVy6XLPrhRxe2a2HHfB6sEgoLKtTNZ8VzPXVtLiVOHbUS+abqwfrE6i9e5uJJ57YkqhxYO7U7FNNNVWuB1NcqBqx51HPd9mx1p616mJw1Yd8SPsK5uWt1RA+SEfRhD5C5tCE9dqEObT7eq3rrJfSRVKvk5odZzYc32655Rb7Lvlss80macO8++67Zvvttzf6apqIXJmVQg6vzMYDVa96sEgoLKtTNZ9VMlfJXCXzcSkQt+8aB+ZM272SxqMp+++\/v5lkkklST6bvvvvOPpPKYywXXHCBWX755WUnWReoVReDqz7kQ9pXMC9v4YbwQTqKJvQRMocmrNcmzKHd12tdZ72ULpJ644SmoTInbpw\/JIzBux0JfbzxxuvU3i+\/\/GJee+01M3ToUDNixAizwQYbWFDPAn\/JoNqlTl0MDjm8JLQLaV8PFgmFZXVC+CDrIcw\/ol36CKFTE9ZrE+bQLmspSVtV11kvpYukXmwGuH\/961\/mmGOOMddff71tY+qppzbzzz+\/\/QPYP\/vsszapDBniKOuvv74hxnzKKaeU9Nll6tTF4JDDS0LckPb1YJFQWFYnhA+yHhTMm7BemzCHdl+vdZ31UrpI6iXmZv\/hhx\/M3XffbSXv9957L7ateeaZx6Z9RbXelNhyf6J1MbjqQz6kfT1YJNtKVieED7IeFMybsF6bMId2X691nfVSukjqZb6ahjp99OjR1smNVK+48pOylXzsSOxZGeIkg2jXOnUxuOpDPqR9PVjKW60hfJCOogl9hMyhCeu1CXNo9\/Va11kvpYukXiaYSxppap26GBxyeEl4EdK+HiwSCsvqhPBB1oNK5k1Yr02YQ7uv17rOeildJPUUzFOoVBeDqz7kQ9rXg0WyrWR1Qvgg60HBvAnrtQlzaPf1WtdZL6WLpJ6CuYK5ZJ101NGDJRe5UisrmMtoGUKnJqzXJsxBxun6Lp8K5lIOddF6dTE45PCSkDqkfT1YJBT+T50PPvjA\/kkqRIaMGTPGTDHFFB1pJZPqzjzzzIY\/eUsIr6V9hfSRRSPGEEKnJqzXJsyhFWsppI+6znrpmCX1VDJXyVyyTlQy\/28OZymxAKlBgwaZxx57TPpJar2lllrKnHTSSbkBPQRopQMv2kfZNGK8UTo1AQibMIeq15K0ferFrVcF8zwUbHFdQutuu+02c84555h33nnHzD777GbXXXc1a6+9tjiMri4GFz0cpSQOaV8PFhmVHZ16rjXQjNdrWtlHCbV+GvWS+W7k8NwPYyQdXEGDifm46Hoqk0ZWgo+hUxPWaxPmIF1zRdeStH0F8zyUaoO6Lovdueeea7bddluz5JJLmkcffdRcdtllFtD32GMP85vf\/CZzpArm45JID5bMZWMrlEmnoq9cdRUwn2zzwWaCWfrLCJtSK45OZfKBrovyIgSkmjAHKXND6BTSR11nvXTMknqNVLOPHDnSgjZZ6TbZZBMbC88jMmS0I7PdpZdeahZbbLFM+tTF4KoXdEj7erBkLhsFcxmJSr3wJAFtE9ZrE+YgXBKxKnDpt9J6qmaXUqrmekjlvK3+yCOP2Jffevfu3TEinmjldbdll13WHHDAAZnSuYK5SuasGdaQe\/dYurzLPICLSoMqmZerIVHJ\/HtTdD9I902IoBHSR11nvXTMknqNk8y\/\/vprs\/vuu5tevXqZ448\/3vTs2bODDrzwxsMxX3zxhX3qdfLJJ0+lUV0MJtseJoGtt95a\/AythNmuTkj7UpD65atPzQ8v3G8mWmCFVJtxUaAKmYOUViF9SOhUNY2YZ8gcqqaThEaMIYROrehDQqcQPjRhDhIa1ble6zrrpXSR1GscmDvpm7fVDzrooE7pZlG1H3fccea+++4zl1xyienTp09bgnnVCyukfenB4kA6yx5aFMxD5iDZGNQJ6UNCp6ppFDqHqukkoZEvDRdZS63oQ0KnqtdSKJ2qnoOk\/TrXawh\/pHOrul7jwPyNN96wTm94raNKjxZU8Hi5A+Z9+\/ZVMB8wINcaa8LhKJ1wyAaX0CkUzCUx2rynQIjcXnvtZcO20krRWPaidJLQKBSkWtGHZD0VpZEPcFIPV4wAACAASURBVNLLjLTe1VdfbQbk2P8hc5DQSMFcSqX4ernAnKdRAUueP\/3mm2+scxm2RBJfTDfddOO8eR42tGJfVwHmkoOw2Gjjv8pzABfpN6R99+0ky21sfpPigYxq9Ns7zjLSenlpHDIHKc1C+pDQKZRGp59+emlx7NAEsN9oo43MTDPNJCWRrVeUThIa0X4InVrRh4RYRWnk01e6l6T1mrbnJHxIWq+OP3kvONI+W1FPBObffvutlWSHDRtmxo4da8e10EILWQezCSec0Oy7776mf\/\/+ZuDAgZlOZVVPqkwwryKpRdXz1\/aVAkoBpYBSID8FiiZmyt9TNV9kgjnS+BFHHGFuvvlmm3iFt8ufe+45w9OogPkEE0xg9ttvP\/PAAw+YCy64wP5eZ5HYzAldi3q6J41Zosqsc77at1JAKaAUUAqEU6ComSm853JayATzyy+\/3MZmE7O96aab2uxp2J2xnwCIU001lVWz7bTTTmaBBRYwRx11lJXW6yrOm33KKae03uyTTDJJx1DyerPXNQftVymgFFAKKAWUAnkokArmSOWo0AFEHxijYO68xJ955hlz\/vnnG4C0rlJmnHldc9B+lQJKAaWAUkApkIcCqWD++eef2yQreDz6nuFRMKfDuH\/LM5Ay6yZlgEPLcPLJJ4szwJU5Jm1LKaAUUAooBZQCVVEgFcydyhpb+WGHHdbh3BYFbqTho48+2j5oIknGUtVkXLsuNzuPrJB4Zemlly6Um73qcWr7SgGlgFJAKaAUKIMCqWDuVNbkND\/rrLMMiVjipPDXX3\/d5kJfccUVRWlSyxh4Vhu8mjZ8+HDrhV\/01bSsPvR3pYBSQCmgFFAKtAMFMh3gXnzxRatqJ2fuLrvsYtZZZx1z4YUXmmeffdbaxwkF453lt99+26qvF1lkkXaYl45BKaAUUAooBZQC3YYCmWCOc9vf\/\/53mxr1s88+iyUM4WmkSV1\/\/fU7pU\/tNlTUiSoFlAJKAaWAUqBGCmSCuRvb6NGjzQ033GDuvvtu89prrxlU8NjSV1ttNbP55psbYvS0KAWUAkoBpYBSQCnQegqIwbz1Q9MelQJKAaWAUkApoBSQUEDBXEIlraMUUAooBZQCSoE2poAIzMnwht0cr\/WXXnopcTqo2skAV2fSmDamtQ5NKaAUUAooBZQClVAgE8xJ27rHHnskOr\/5o3KPr5DiVUsYBT799FN7MSKFrgsJDGtRvy5CgY8++sjgBJr3JbEifek3SgGlgFKgKAVSwdzlMv\/b3\/5mH1NZa621OuU6j3bao0cP06tXr7Z4CrUoQdrhO4D8kEMOMQ8++KBGCdTIEIB8n332MT\/++KM55ZRTzBxzzFHjaIp3zaNIb731lvn444\/N9NNPb+acc85S92jV7TPzqvuouv1WzKH4CmmvL1vBi\/aacTmjEaVzXXDBBTtlgCuna20ljgJ1ATnRCbfccovNlHfkkUeaSSedtDCDmAPZ92688Ub77v0SSyxhtTvLLrtspaGLZc6BydMe7xvz0BBP\/JYN6E888YTZbbfdzLTTTmtfG1x11VXtY0W\/+c1vCtM++uH7779vEzmhYXMFTQ9pjaebbrrgfqpunwFW3UfV7bdiDq1YS8GLRdBA1bxoCp3iSClK58pB5udmF\/BEqxSgQF1AzlDvvPNOs\/fee5v55pvPnHbaaTbssEh57LHHrBbnk08+sRIgBamQZEOk\/PVfsSvSfto3Zc3B76NqQP\/iiy8MNCNb4cMPP2w1WzvuuKP505\/+ZCabbLIgEuHjsvvuu1sfFi4Ns8wyi71gkeypDH5U3T6Tr7qPqttvxRzcIqlyLSUtRKTor776yvTs2TP4tcxW8IJ51EGnoI0s\/DgVzLEVkmv9hRdesFJJ6OEiHFO3rNYKIB87dqyVjH3JDx5\/+eWX5uKLL7bS27Bhw8wMM8xQiAfOvwJHSPL3zz333LYdchT89re\/NZNPPnmhdrM+KnMOcX1VDeiuzw8++MBqNAB2QH3IkCFm5ZVXLqQSJ50xPhdcpM4444wOKZy5nHrqqeaKK64IenCo6vahSdV9VN1+K+aQtDfKXEtxfQDid9xxh31N85\/\/\/KchcRjat5133tn+d97SCl7EjalqOuWlQ0j9TAc4JKw999zTSlkwq+hBHzLIpn8rAXKkNzbP\/vvvXwgU4SPfkm4XPjpAB4D5d\/jbr1+\/whoYbrt77bWXfdseQHJA3grelTWHtLG2CtC5mJAqGU3Yu+++a3bYYQfLr7waDd4j2GabbaxEjhOlX3iqeKuttjKHHnroOL9J+VV1+4yj6j6qbr8Vc0jjV1lrKdoH7ZJA7PDDD7d+VJjPHnnkEbPeeusVdtZtBS+SaFUVnaR7qax6mWDOIXbddddZuyG3p7Si3uz52QIIkir3nnvuMQceeKBVsSI9+8WprpFu0ZTkBUoH5M8995zl4x\/\/+MeOPniz\/ogjjjA333yztW3TPjbcvAVnPcbO99h+o4W1gxr5qquuspoA7LYASp6+cMjkz4QTTtjJpl\/WHLLmXAWgc5DwZ7zxxuvUPeuCw5ILHGDuX8Cyxsnv2Aa32GILM3ToUHvg+gW75LbbbmtWWWWVwpe3qtvXOUi43LlOVWspOpJXXnnFXjIxy2200Ual+MG0Yj1lUTR0z2W1X\/XvmWB+0003WZDhIMsqCuZZFBr3dzbgQw89ZLUfFP91Ov6\/SiB3o\/HBsAhw0A4XPlT0vFLXt2\/fjomijgPoUfm+9957nQiAlgAVcFbYF5oLLgnXXHONXYdoFUgjDNi5y0AZc5Bwr0xAdxIOqr44sHZzuu2226zKPc\/B6QAbHwjeTfAdGj\/\/\/HP7eNKAAQMKg3nV7cOLqvuouv1WzMGt2SrXUnRfsBfvvfde++BW3GX822+\/tcIBe1R6WW8FLyT7O2TPSdqvsk4qmDOxfffd19rdsI0stthihex3VU6gCW0nAXorgDwK6ADHwIED7Z88XtVIkIRxnXvuufYpXMqoUaOsJmDEiBG2rY033tje6PGiRk0HQPHePCrlpL6cUwwXgdVXX906cT3++ONWGvBtwfTnb8Qic5CupbIA3fmkcBGKXoLcWJyZC2koz6uEjJGLITxAne4ncsLGCZhzOdhuu+06TdtdlrJoUXX79F91H1W334o5+GAOyFaxlqJrAX8YTFsXXXSRicsp4rR0OLxGTTyurWj+hlbwgr4RLvAdmmiiiRKXeNE9l7Vnqv5dFJoWcoOvegJNaT8K6Kigud1WoVp3NCN+GpV1GYCOfXennXay0veGG25ouJ2Tn4CNs8IKK1hVPkDsisthwKbGuzoua6Br8+eff7ahVO553bQN2dUA3c0Rs8PBBx8ce6nhQrPrrrvayzQHZNpB5O8H6EYZf\/zxO20Tni1GzY42wD9suTAQyQD49+nTJ3NrVd0+A6i6j6rbb8UcHKOqXEv+YsBZFvMNl0vWZLSgruYyjaY2LgoqKX9DlbxwYatcRHj9E8dSzqQkrWDRPZe5aSqsIJLMcY7S0LQKufDfpn1AJz57rrnmKt1G7mYxcuRI+w794MGDO71BHwKG2OQHDRpk3nzzTdvNPPPMY\/7yl7\/YGOqoTZiLBGpytD5xYM7mY+Ndf\/3145gesjgRMoestv3ffQl96aWXtgdc3uyH8Pzyyy83Z599tpV05p9\/\/nGGQB3ogJYj6QDNM+4XX3zRgjmXB54tpgDkhLFRkPBQzxctVbfPuKruo+r2q5hDq9aSuwwStRLVjjEvp9FNwo28mq1QXjip\/\/TTT7chtwgwnDtpuSPK3nNF91Ke7zJt5hw0qG44cLNsm3k61rrxFMiyoWfRDYkYkESq32STTWwCmKgkB5CjEmczxiVCCQFDNg4XEUpaNkCkCFS9SKSHHXbYOBKps6HhTMcFIeoUiGTPc7znnXeefZK3d+\/eNr6d+GlU9iFzyKJxFNBvv\/12e3HJAkB4C1+QTABTYnMp+AQApFze4F2c5O2knYUXXjiWHnnGjIqUZ4tJiIPWrUwgZxxVt9+KPtp9DnWuJWefJ0sldvHoexyo2fEBAui5yMeVPIAeygvnsMeYMPVxlmDz51whI2KSU3GZey7P\/ixaNxPMORSRCJ5\/\/nmr6lt00UUTw2Q0nWtRNnT+LhTQnc0H2zLqLoDCAUQWkLuROL7j+V52bngHtNjSk+zASargOIc6VMLEskM3VMRbbrml3bBVzqEIp52EgZpv1llntcDttBb4KjD2NLs4F2uy9HGxnmaaaYoMwX7jH45cuMqSyN2Aqm5f5\/A\/zURda8lJu\/hlIOThC8OaJqzysssus5fFrAgMKaAXXU8u\/8T9999vo4Uw1TlHUP+MTQP0svZc4c2a40ORzRz1qaSoN7uESrI6VQA64R9pErlsZPJazAEJe8YZZ+yQvJFCcaYEvNjsXBDjnN+ckxaqOhzlSDjjO9SRmILQK+z05D4gvh3JfMyYMdaZTGLzlc+knJocXscee6ydOzQhGRO2O0CchEzE6TOXJLs42gzixqGH8x8oMjJ3ONLfrbfeapsIVa3746i6fR\/Mu+sc2mEtRRPHwBeXPAatmyQ3ggTQi64nd97hUR\/n9yUB9LL2XJF9mvebVDB3cXeEzkiKPoEqoZK8TpmAvuSSS1o7UZJqXT4qeU20AGgGAK4\/\/OEP9qEPntJFRY6KC1Vzkhe7b+eaeuqprSTK+Pn33\/3ud9bWH00562Kynfo4a6SOvkjDhG+1IiESc+ACguMZTn8nnnii+frrr+3FBgdBTA7R8EQ3D+dnsNRSS3XYurPmGPc7zj3QngsTNCwTyOmv6vZb0UdXmEM7rCV4wZ4EI0jrSg4MCYj76zIL0Ivywmnm8DcBzNlXUUfbLEAva88V2ad5v8lUs+dtUOuXSwF\/sSG94UmaJ2mMr3Ln1S8WNPbdVhQug3hHu\/hw+qRv8hb8\/ve\/z0w2QaKZK6+80tre3IMtqITJOBV1qKNtl9wI1ZhEciUTGlnSUL0xTpLmFCnQmL4Js\/NT1sI7vPux5zszh3NSInkO6nKkBjQJxOhPPPHE5vvvv7cheGgvoocikQHY52mPw6lo4ZAj6x+mjLKBnDFV3X4r+mi3OZS9lqAhWhnCREPWUtE1GP0uDdBDeOEDelIODXfG4gPAvouaFduJTmn0zvRmf+mllwyvpuW9bZXFZG3HWFswr5kBYEU2XpoNvRX0dZnbkMK5kESd2coYA2Dv1OwSmzLS8F133WVBFC\/yohccl10PUI6+NucSLnGxwHOcfcTcncMb3rQuHM03ISAtk5DD1zwkpeMtSjsuMjjg5bkY5umr6vYZS9V9xLUP78hk6Bwt89Akrq50DmWupdAxV\/V9GqBL6RQ3NimgI4VLwz6rokFIuyKb+dtvv23WWGMNK3lw6MVJRSGD0G+rp0DdgF7lDDkEAGVywnOz3mCDDVK7c+pJTA6AcNHX3NKAnAGgWcA+TsgaKm1nH0edzgtvxLly8fDj53FaZEzuxTnaSUvHWyVdGT9hPGVdvpLS8XaVOUAPeAawVpmUKI4eZa2luLZdDHYZzx+H8jJL5Z7UPt9xOSfEk+gWBAf2G1pA9psE0EPHXvf3qWDOAsIOiQoRZwIK0gL2PuJTpan66p5kO\/Rf9cthkjk2AdBxGuP2TMpYJ+U6VT6XTaTzrLfYHfhzAcCen5R9LY2mWUDuf0u4IH1ccMEFVo2OExvZ17CXA3BpLxLWAeTQB49kF5I622yzSZZXYh1JOt6gDmI+LnsOrgvnaPn000+3HNAZQ8haSqJxWU8HJ5mb8vI2b\/4G\/zJPqBzmRCJH0KCgUSafBu8QQDsis7ChF01bnXcurawvspmj5sV5iZsPhHA3H7L\/4EiDQ5Kq4ePZ5uJBeUwFmwyZ3eosRQGdfN5kaKrzAucyXJGUxk\/+AD3zvi7mAzoX06iKvAiQ493L3sDJD\/s2cefsDT+eHFBHrY9tnbA04l1xvltzzTXH6bJqIHfSClo3Px7YpePEpINDnp9rP8\/adbZOLmBEGWSl483TtqvbijkAAD6Nygb0uDmUvZbiaOtCt8p4\/jjP5VbCZ\/anNH+D03IReuZ8cZgbFy6SVrH2nFNpkyV0EZj7xIdIJL3AKYD82jxdhx1U1fDxS9QlP8GLkuQK8847byG1JYsbZyU8rkkGE2LqcICOuhm1VJbjl4tN54nDtDzqWZs07j31rG+ia498B0iyhKtQuFAi6SY5xaW17wM6NCVELEuqTzq0OCQYFxKtX4i9JZMdNHbqav8NZcZApAEmAt\/TVgLk0XS8eWiJtEycO9IKhyCOUOztV1991eaUwIkwJFGUywLG05hcADE1SNLx5pmDy1fAXKqaA2p1v\/2yJfS4OZS9lpJoWtbTwWUDed41wPshOIXGJZ9yaVmJhnEvQkJftGIkpCo7h0aesZddNzeYRw9XpCVuPdiRKNjUebwBScNJJGUPuqu0h9MG3swc8ABFXB5jyVwc6KAZKfIEalwfHFDcWLOyluEohgc50gjS6zLLLFPoIkE7hJMh5RGaFc0XLqFDFXXyAHrSoeXauPHGG+0hgUSLZIXUjc0OuzMSuZ+qFeBk7wCYaGt8PkiAPCkdr4RGaNm4vLAGAO3FF1\/cfkamLMbChYZMc0ikRfYwhyUH65NPPmnDlbgwln1oxoFg2XNIAvKyAD1uDmWvpbT14PoPef44DcjL9rmImwv9ox1ed9117bkSVzAXoRnFrBZ9DliyX7pKnUJgziLAUYfwGuwSLECSBaBG4zf+LS3vbVchTsg4kXjwVAaAkUxIOxqXdzurjyqAPKtP9zuOjyQMwm+CzbL22mtLP+1UzwE5a4OEJ3Wq6uMmIAH0tEMLECTEjQuboxFgPnz4cHuBQSODjRzNhP+wTdxYykjHm8YkB+RoB5Bodt55545Yf98kRJhcnlfaXJ8OINBIQTMu9UXS8abNwWmWcCpEW8TDPiQJQvNR5hziJPLouIqq3NMuI2WtJclm9QE9rx057XJbps9F3DycDxJ7lwQ1CAlJyZaSsklK6NOV6ojBPA7A8RjEJkgmLh6acHZzJFIWBlJcnleeuhLhJGPlMMBOTtpSNmjSq1hJbSUBufNsnWKKKazHZhUF4CG+m0sIUjTaF\/e0aZ7+2h3I3VzSAB1gxs6NPRxwjqrikWxQsTtHuiiQAzb33XefvdhJbPNlpOON41EakLv6Phi69+ZRwUuKD1B4e+OMGH2ZLU863mifLh4YTRGXHlc4h5D+3V4oaw5OdZ8197yAngTk9FP2WsoaO7+78eR5\/jjtclumz0XUXyHqg0Ta2KwHmRyYsx7dw0ISunS1OqlgzgGHXYWMWs6phwkidaOOY\/PEPV3pnrckdjbpzduuRqii43WA\/sADD9iLDYn+JfbuJCD3M6ORN52FXNT5ED7h2Jb0gI5vu4PXaV7XcfRJAvJ2CoXxx+38EqAHkrQfksUhwp+kZDUc\/IB5v379xpHI6QOPWvxMpJ7zcYAeko5XAuQhgB4FKJwlkZjKSscL7QEbkt2w3tHwICy8\/PLLNtJmxx137HTJKgLoaSCbtf+lgJ7VByrhstdS1tjzAnoSkJftcxFn5ojzQQKsSQtN8Z9K5v8ZE07bqNile09Cr3asI4ozR9WKg8Gf\/\/xnmyyBcJ602FNuzdgKUb2jYsxyKmpHwpQ5pryAngTkLEzs3NiGCL+QSHlJ8wDIcUpC2iS0K+5S5ja5c+yKPtqSRqM0ibysUBjohDOmU7GWybM8baGJ4nKLhorYdYDGqdbdPuGQznug+IAeko5XAuRk60O17i6GecAwyfaLNodnJ8tIxwuNOX84UxAuuDRlldA5ZLUf\/T0L0LOAnPaqWkuSuUgk9CxzU1k+F3FAnuaD5FJHY8YiXz\/ObRQ84jnncPaMXtAlNOlKdTJzs3MIEbaiyWLC2CoF9DQgJ3oApypilQH0ok5FDsjZHKeeemrmoyS+hC4B9CQgLzMUxsWiXnvttVbFGhoLnYe7qIoxQ7hsUdAThy+cQFH5Aub+wcFYAXcc3nBg9FO+ZvXrA3qRdLwSIHcev5hRkHz91LOofXEeSlK5pwFUWel4v\/zyS+uwB1ii6YvLguiibIiuwWmPqBEuJhJAl4BsFp\/c7w7QiRQ599xzO5xek\/po5VqSzCEN0LPMTVX6K0h8kBA68c8gdNUVLn+Ep5GDIukdCAldukIdsc08azJIeBzipKzsyinxsuYZ8nsWoKc5u\/k52kldys2VW3De7Fx5gdzNVwroaRJ5WaEwDshxsomq1UL4k\/WtMw9g2uDpSUwP8AG1vANENCf4RuCYxT4A0LgQE0te1O+gaG4A5gP48VIbdm8uPdGc9W7c1I2jpX9AA6JRLQ5OkqSe3XvvvW0fRUtaOl4unajRuWwwF\/9QZnyknEb64qlLV3gbHq0AUQKtmoPrG+csNB1+hEKUTnWtJQl\/2OvYquOeP04zN9G25PKUNIasS5XEB4n9BqjjoI0GGZ+udnO4lfCgSB2Rmh2JAE\/dJFUsHXPAcVh3dxt5FhPcguSAisZ48xuHEs5ycTmzfUAnbpLY5Dy5tYsCuRTQs5zdygiFqRPIndrYT1jjR21wYONkw2FC7gW0BZgBoAuHI45wEn+JuDVUFND9wxUTgA\/YWUDuxkEbqCtRl2OnLrv4Gbzi0vHi14H9HUHBjyWOZpYjlh8zIHZV7KOEggLoHOZVzyEPTXy\/lzrWUp6xFqlbBNDTvPs56xwgZwlERcbblG9EYM7hxFvlAHYSeCiYy5cECxKnDSIA8krWPqBzycKezd9ZJRTIo4BOWCJ9++Fq9IGTFxJU0m04JBSmLiBn7oSfkWVuzz33tE6M8A2TB2q96aefviP+H3Ugqv9bbrnFhmZh6yb8C9DPy+soT8sEdFTRSU5DWWsp5Pci6XjdM5REA6D1QNoiOgCgRkMCiO+zzz6W1lyWnNMT5gG0BnmjMKpOvdwOaymEh5Jv8wB6EpBz2eQSx4XHjyxQQI\/ngAjM+RSpg4LTFbffqIShYC5Z4uXUyQvoZQG5D+hPPfWUzbpWxA4lcbSJUqoVQI6dE2D2HTadnR8V7j333GMPFfe7zwcf0MvhcnwrZQC6i3VHLd5KM0VIOl4OdkLSfHsomgLU+zysE43ocLH\/xNLzloSkRMOeQlIvd4W1JKFJSB0JoKep1l2IJ9oYTDwK6OncEIE5mwH7D1IIGzIuD7aCeciyz\/+tFNDLBvL8I43\/Ig+gtwLIUYXj\/IXk7Sf3ceFgaBo4UEhn65c6AV2ajteN1z9cuYyTra5oIqAi64D+Q9Lxom5HOudCg7S20korWXNGXCGMDZMHGe6kzwaXlXqZZDY4QOKkiuOVK+24lorwMc83of4KCuhyaovA3B1iLHYOM2zjOP\/goUuucIqCuZzoZdXMAnTUkyR+eeyxx0Re62WNS9qOBNBbAeSM16lycRDzJTnnDESsKvsAu3nUd6QOQJem443yIs2GLuVbu9dz2gvm6vJxZ425rNTLrBf8izgjURHjuOdKu66lLNqE\/h7qryAB9DgfpNBxd7Xvc4E5k0OCIZyJGy8eoyTDwJ7OpiGxjDrAtXYJOCDhMoUNO5prHU0KalU8O9uxZIXCAKI8HdoKdTCvR+EF66vSoZl\/CCelvHR8IOMfTlxFwwZbwaMmAzpe5MQU33333TY0TGIvl4Q9Sfjiq4zJd09Cm6avJQldpHXQXOCvQP70qBk3C9CL+iBJx9YV6uUGcyaFtES4DZI5iSawoxP7zL8pmLee7RzOSJZdNSQwLRQGCZQLJJnEqi68gbzLLrtYx6roozhSQJfwoRXpeLNo1dUBnfGTdwGtE9Ef2P+dBopY7zzZFqGVJOwpjaZR2y\/agVaspSw+d5Xf2ePEg+PYmMQ7P+Uyl2Vput2uQoPQcRYCczplM7F5sDOycMnCRTIMBfNQluj3dVGAA9klTUG6jjpVSQA9a+ytTMebNZaigN4O6XiJS0eY4KU6v+AUR0w\/9vS80QNFvaTjnLhasZay+NvVfnfZ5TBRJAE6Ejg+W9TF7IVmmNA1LcYUBnNHPN+OjrpdwVyXVVemgHOcIjkKL35FSwigA56tTsebxQsH6AAg6mkShWSBYFnpeLPeBsgaO6pXwmZ5M52MfJiYiLBIcorLas+X0KVvKaR5Y1e5liRz6Yp1sgAdvMEngVcciSpZZZVVMtdrV6RDkTGngjnqQLxPkbhJ55pUUJFgM8eLk9tyWnKZIoPUb5QCraIANldyO6OqTUrKUwTQnVq4jnS8WbRjbKR9JbwuDchdmB6+BUhP5Jp3DrBZfUR\/l74NkLfdMupLJfSsjGVVraUy5lh3G1zEyBrqHvDiIkb+APIf+IAOcPMyGiGwLn8A5pXzzz9fcSbCxNLSuda9OLR\/pUBZFACoSPTCAUM0QFyaUg5yQrt40EHi7OZ7vNeRjrcM2pSVjrddwyV9GmUBehaQu7aqWEtl8LLONvzU0P44Zp11VhsVRRIg6I\/JC\/rxwM4aa6xhSDiEPwu+CPzJ0iDVOcc6+s4F5tjKOIg4mOIKxO3Vq1fhlJV1EED7VApEKeCkaDzTF110URslkPRMbB7q+YBeVjreVtqvy0jHmwTkaAGJuminAzot9bI0J31VaynPumunus6JDV8HLsPYvZHS0Vjx0iZrgFS85HrwnRLdHIjbJ39Ad3+JMxZ\/f01C5v\/Whvik6ORAQ42eVtRm3k7bRscSQgEX305OdXKsIzEsvPDCwWCTlRsgacxJIFiW\/RrVJtEoHJZpr7qFpOONmwN05sEc+kZ12srX7yTrIyT1smu\/qrUkGX+71XGZ+VCfu4RFLuxs8ODBHU8H4wMBsHPJI2SNDID4Q\/Tt21eFxQSmZkrm7rBgQSKdpHkO8pCD2szbbfvoeOIkb2KLs273AC82PXJ8f\/XVV+ZPf\/qTTSka+gpTXkCPA8Ey7ddpb1THrR5Jsp\/od0mXkQcffNAmVkFCI20nh3WR4vwYFJ6WgAAAH9lJREFUUMfyZHMVhTlcc8011skXwQZ1MHkHAJksjUJVa6mKeYa2CS+4BHExi9KFJ3URDJG+4XUUyHmMiPA0Mv0R8py1R0PH2qTvU8GcxYvqgzzcvOZFIoSsRdsk4uhcmkcBDlVUegAYtjhJ8V\/nIlwNhxyAnXDMovtBCuhJIFiW\/TovkDt65QH0pMvIq6++ah1sSUCFVF7UlCG1X0t4nVTHJaMhPz+XBUCGix6FrIASvwnqVrGWQuZV9rdI0jwLDBjjb\/L73\/++0x5BA0N8OGDer18\/M3z4cONL5IyHcDO0wQ7wyx5jU9sThaZxA0U6KXpwNZV4Oq+uR4Gk5zQlM+Eg5ilQstLhiUt+cFK\/opoukrAnC9DT0vGWYb8uCuR5AD3pMoK6lYdMAEXMF5gz8GTOW1oB5GgleUcdLQJROy7LontCliiAuDS\/aXMpey3lpVuV9eEtGizCMKOATtpV+L7FFlsYNLlDhgzpUK07fAHwiZRQMM\/HpVQw5zbKYwFI5NEHJvJ1o7WVAu1BAadtQr1HYpgiAMJMCMfkoZMpppjCgnrREpKON8R+HQrkEkBPu4z4CWvIInnppZca8uLnKa0AcsYTZ+d14yRMj3XEk8AIPUVKWWupSN9VfZME6Ow\/zCk33XSTDTcDzHmUxgE5FydMtaSh5uKU5r9R1di7arupYA5hcfx57733rJ0jJBlDVyWQjrt5FGAtP\/TQQ\/at61D7dxnUCUnHm0fd7cZaFpD7gI5kTcKZqLo57W0AH9AB8qQwwDgatwrI6fuOO+6wAMOFI\/r2AeYOzC5IkuQU1\/I\/CiQButNoILnzPv16661nNVuo6JHKSWCEpkOSV1\/p\/T8KZDrAOcJzeyIHcjTFpRJTKdAuFODy+dFHH3WyZQMYvIGNBE2qTx5wiDrhtMv4i44jD6CXDeRFx+y+KwLorQRyxjly5EirGibGOfpkLOlFt912W+sIJ303PZRmXen7JEAntI8QMzL4ISTiLMfeRUvBxRBHuOhjK11p3nWMNRXM2TTYB59++mmDV\/sEE0xgX0ojZCCuqDd7HSzUPqGAi19FekLiRuXpnN0GDRpk8yOwfpdaaikL6n\/9619tnOvqq6\/eCI9ZCaC3G5AXAfRWAzlj5CUvkghhDog+qerAHGBaf\/31dTPGUCAJ0Ikoufbaa80tt9xiHVKXXHJJS+f+\/furf1aBlSRygOP2JCkaZy6hktYpmwL+a0o8wrDPPvtYDRJ2cQ5inqJkDfO8KR7UALtfkAzmnXdee5iwhklYkZXatOw5lNFeGqBDC0KqcNzDe7jdQn4kEnodQO74AiCRpjqavtap2a+++mobXucK4MQ3UW\/uMvjcDm3AC9YR62nppZe2ewcTBPsGW3i0pDnFtcN8mjCGVDDnACC+NiOvTAcdNANcE5ZE15pDEpCnzYK863jbEgpFhje8lJ999lkL8kjthGFKQ43ajVou3jrOfs0+5k+7qi99QAcYeezGf+dBmnWtlTyJA3OnAQHAzj777MKOca2cR96+kKh5mIf95xeAnHc84gBeAT0vlfPVz7SZ52tOaysFWkeBIkDO6HC0wfEG8D711FOtlOqSsPBb0cdDWjfz5vYEHwj\/41K1zDLLtP1EcdjCcY8wqrnnnrvjkRA0Qc4psCkhvfCGPYe5ir+ZMzHhU001ldltt93MBx98YB599FErrfsg7wM8GjNyCmCqjYtDb3uGt\/EAO4G5k8THH39865QgXYSEoBDkj\/PCJpts0nYqvDamvw6tIAUkQM56pMSFt+DRTtIPVM+9e\/cuOAr9rLtTwI+J5l0KHgdpIpDD5\/vvv9\/gI0C4MsUHdDQpJIPhUSL+nVf4kMSdeQvNFyla\/TLjjDMaQvu4BGkJp0AnMHcJNWaZZRYbO+l7rjtnOA7GddZZp5MTnIvdHTVqlL5nHs4TbSGDAhIgR9W55557Wi92Ql2iz\/ISbkR2Q7KP5Y1vVgYpBRwFHJhzXiJxNhXImS8mBcDXabPSAD26QhAUP\/vsM+uzQkZRnKrx\/l9zzTXFQqOuunQKiME8DegVzHWZtZICPLyAWu\/777+3ar7VVlut04HggBzJgN9XWWWVcQ4Ml4mK+GH1Qm4l95rVF5dCwtJIPoQvQtNU6z633nnnHbP33ntbQc+Pt0+S0JvF6fafjYJ5+\/NIRxhDAWJ\/8VrHxOMDtgTIae799983e+21l9l33327rLObLoz6KeCyw3GxbDKQQ2m0s5gRCB3DgdQ3wyqg178WFczr54GOoCAFooC+4IILWskhTSJ3XXEwYdebc845C\/aunykF\/uNMiaRKVAQRBFI\/o65KO8xSV155ZexztQro9XJVwbxe+mvvgRRwgI6zG0k9KEmq9cCu9HOlQCwFsAcD4k0HciZPet7tt9\/emrbQakVjyhXQ69skCub10V57LokCDtBxsCHbG1mk2jWWuqQpazMlUoC3yQmFS8psWWJXXb4pwJq0tkjoSU+\/KqDXw2YF83rorr2WTIEkG3rJ3WhzDaPAY489Zq644gr7xGnIC104AeMkXPRN9q5EVpIuEZ6Gn8A555wTG1rmAB0Pf5xM11hjja40xS45VgXzLsk2HXQcBRTQdV3kpQD5Bu666y6bAKVPnz55P7f13ZvtJEuJZq0r1GAX+Ii9NnDgQDPHHHPYFzX5O1oAdLL29e3bt1uYIOpmm4J53RzQ\/kulgAJ6qeRsfGO8oIdUDpiTkz9vcUBOqCPx10UvBHn7rbu+n3p31llntX4qvGugpT4KxII5tymepPPtjoQA8ZQkTkbR33AAwTGC2xkZtUjvp0UpUBcFFNDronzX69flG8AOnPf97CYDuSSPP3XIonjQQQfZ7J\/E2\/McLDH3WlpPgVgwl76SFh2uvprWegZqj\/EUANCJicXzlgMm7iUnpV33owBpRckxsMACCxjSiRKeyPrgtb3NNttMTJCmA\/kNN9xgc60D0Fl7h8yfxNiPGDHCIKUfcMABZvnll++UQVRMWK1YmAKxudmlr6RFe9VX0wrzQT+sgAKffvqpTeOadRhV0LU22YYUwNmNd8c\/+ugjOzrWBWpxQIt37nfccUf7FG7WmmkykEMXzn\/ebSdVLeYHbN5ZBe0smRnRzD7wwAMWyDfccEOz1lpr2YuTRgpkUTD8d301LZyG2oJSQCnQJhT49ttvzY033mid2nDAIpXvRhttZE1\/7mU8TIIvvPCCzaOOpE6a0miZffbZzcILL2yWXXZZs+qqq3Y8HtV0IHd0gEZ4rPMUMC8M5rkQ8xohwM7TwtD20EMP1cdUWrA\/FMxbQGTtIh8FOIRvueUW+5zikUceqa\/w5SNft60NOJNnAN8ev2y55ZYWUHi6M64ghbLeUBWjgudBkccff9zagf0Urd0FyJ10Tiw577EjbRdxDuy2C7GmiSuY10R47TaZAnfeeadNy8pjDoT6ICUVKeRpR1W49dZbB8UQF+lbv2ktBZwK\/be\/\/a0ZMmSIWWyxxWyq1UceecR6WU877bSJA0p6QY8nO5FIMR\/iAMz72\/TTNK9155mO+cF3YMNMRQ72ueaayxx++OFmookmai1TtbdcFFAwz0UurVwlBZwalGcWkY6GDRtmZphhhkJdAuQ4wHGIq3RfiIRd5iNAh0dzcGZLSmKSNhnn0U6OdWy8SQXVM7ZfHOeaVF588UUL4mRQxIEN4MaBjWim2267zWo1Lr30Un0quM2ZrmDe5gzqTsMDwAFgHj\/p16+f9YotUhTIi1Ct634D4GDXRV0OCEULzlkkdCGMigvjkksuaZZYYomO0FvSuRL1sNJKK1kHue5WMGsRaw8duajgT7DyyitbEJ9sssnsRYlL9dFHH63SeRsvDgXzNmZOdxsaL5kdccQRhkQeHLYczmnq0Tj6KJB3t1VjzAknnGD+7\/\/+z1x44YWd1gvA\/fzzz9vfuSj6hfz9mHJQHbPueDSE\/0Y6745x0m+99ZZ1eNt0002tN\/+JJ57YETs+yyyzmMMOOywxF3v3W3HtOWMF8\/bkS7cdlQ\/oxLhK4lwdsRTIu+eyIZ0ofhY4apHQioLqncvgNddcY73aV1hhBbPNNttY3wls3w8\/\/LD1x1hzzTXt74AXIVRNe8aU\/QQtoIt71Y1LznvvvWd69+7dIWnzbzi8XXXVVfZ5Uy7RhKVh6iJRGHnYV199dXvZIexMS\/tRQMG8\/XjS7UfkAB21H\/mf+ZMVGqNA3n2XDSFQxIgvssgiFrCRMi+77DJrA0bDc+CBB1onOJfRkodCSBLDv\/EICHZwAG+aaaZpVA5xHADRdBGmxwXm97\/\/vZ3fTTfdZGkCvTBPLLjggvbfncNb\/\/79O8LR\/IQwOKKi\/SjqkNp9V2hrZq5g3ho6ay85KZAH0BXIcxK3YdX9Jzf5bwpPmuJ\/sc4664wjSTq1OhJm0yXNV155xXqkjxkzpgPQ8cznojx06FCDv4Czj6NOR8PBBQDpHLCn4HNAqB7SOv4sWtqTAgrm7ckXHZUx1pbJwZImoSuQ61KBAqiJCa0iSQk2bzK5JamDATgkeOLPAbqmlzhARxInwQ6q9AsuuMCq0XfbbTebYAeTAzH1mC9wgNPSNSigYN41+NRtR5kG6EgM2EnxVB48eLAml+m2q0Q+cbeeyN3fnZKhJAE6lEO9DqhDD3wKiAi49957zXHHHWd9CrR0DQoomHcNPnXrUXIAu0xcpJf0i+R1p25NvG46edYF6mQ\/0QkZ3Uj4gqMXmeJ22WWXRtnIs1idBuh8S4564vSHDx9unQIJ4eP\/8W7X0v4UUDBvfx7pCJUCSoGcFEDyJk8BNvMBAwbYFK9XXnmlwflN6lSZs8suUT0L0LkEkRwHm\/lWW21lszBq6RoUUDDvGnzSUSoFlAJCCgBIPPQBmBOC5QpOcYcccogFeOfZLmyyrauhucLMhLlp6aWXthI1IDz99NPHRoFkAXpbT1YHl0gBBXNdHEoBpUAjKYCaHeDiFS8cuQC4Jj7Fee2119oUrM6T3zGTcM555pknFuAV0Ju35BXMm8dTnZFSQCnQjSjgh+bx1Cte6di\/eXUQad0HeR\/g8fZHnc4Fx49D70aka9RUFcwbxU6djFJAKdAdKeADOj4CJ598spluuukskPMADZI4T8QSL84b7rwI5xdysvPA0dxzz90dydeIOSuYN4KNOgmlQHMpMHr0aMPTpoRNaUmmQBKgR78gpJPseK+++qp56qmnzNNPP21zshOG5lK+Kp27HgUUzLsez3TESoFuQQH3zvZBBx1kHdfwri5akEQBqqy0wEXbb5fvpIDeLuPVcZRHAQXz8mipLSkFlAIlUuD999+372wT50wOdbK6FZEciS\/H27tv3772ZbDxxx+\/xFG2X1MK6O3Hk1aMSMG8FVTWPpQCSoFcFHjmmWdsWBkPpvCu9mKLLZbre1fZATnhW0OGDMn9pG6hTtvgIwX0NmBCi4egYN5igmt3SgGlQDoFyBnOa1689sXjHuedd56Zf\/75c5OtuwK5I5QCeu4l06U\/UDDv0uzTwSsFmkkB8oVjJx8xYoR9FAVwz2Pv7u5AroDezH2RNisF8+7Hc52xUqBLUMAB+gMPPGCOPvpos\/HGG4sytymQd2avk9CJKcf3YI011ugS\/NdB5qOAgnk+emltpYBSoIUUyAvoCuTxzAHQ3377besEWMSJsIUs164KUkDBvCDh9DOlgFKgNRSQAroCeWv4ob20JwUUzNuTLzoqpYBSwKNAFqA3GciRqm+55RabnvXII480k046qa4NpcA4FFAw10WhFFAKdAkKOEAnbO3ss882SyyxRMe4v\/vuO3PrrbeaFVdcsXHhZ3feeafZe++97UMxp512mpl99tkL8Qsaff7552ammWYq9L1+1N4UUDBvb\/7o6JQCSoGIhP7GG2\/Yl8CabvslA96XX35pc6b\/4x\/\/MMOGDTMzzDBDofUAkA8dOtQ+vMKFgEQ8WppFAQXzZvFTZ6MUUAo0hAIA+P7772\/mnHNO069fP\/s+e5HigByNxqmnnmr69OlTpBn9ps0poGDe5gzS4SkFlALdkwJkrTviiCPMzTffbE0KZ555Zm4TggJ591k7Cubdh9c6U6WAUqCLUcAH9D322MPwR5o8R4G8izE7cLgK5oEE1M+VAkoBpUCVFHCAftttt5mBAwfaP1mArkBeJUfas20F8\/bki45KKaAUUAp0UCAPoCuQd8+Fo2DePfmus1YKKAW6GAUkgK5A3sWYWuJwFcxLJKY2pRRQCigFqqRAGqD\/+OOP5owzzjCPPfaYeq1XyYQ2bVvBvE0Zo8NSCigFlAJxFADQjznmGPPHP\/7RLLfccp2qvPvuu2bCCSc0M844oxKvm1FAwbybMVynqxRQCigFlALNo4CCefN4qjNSCigFlAJKgW5GAQXzbsZwna5SQCmgFFAKNI8CCubN46nOSCmgFFAKKAW6GQUUzLsZw3W6SgGlgFJAKdA8CiiYN4+nOiOlgFJAKaAU6GYUUDDvZgyvcrokrDjwwAMNaSclZeeddy78EpSk\/SbW4TnMBx54wD66Mfnkk8dOkTert99+e\/Pcc88lkmCyySYzCy+8sNloo43MqquuaiaaaKK2JdcPP\/xgbrzxRvvsqf+W9wknnGDX2iWXXGL69u1b6\/gdzRnERRddZKaaaqrU8bj6s8wyizn++OPNJJNMknv8ZbSRu1P9oG0poGDetqzpegPzwZxDNwls3MxWW201s8suu3S9idY0YuKL9913XzPXXHOZQYMGJb7n7Q75N99806y\/\/vqxb1e\/8sor5sEHHzRjx461j3fkecCj1dM\/5ZRT7MthUdBWMP\/PpS3kQtBqXmp\/1VFAwbw62na7ln0wv\/rqq82AAQO6HQ2qnPAbb7xhtt12W3PUUUeZFVdcMbErB+affvppqtQKoO++++7mk08+MZdeeqlZbLHFqhx+4bbbCbSTJpFXMi9MDO9DlczLoGJz2lAwbw4va5+Jgnm1LLjjjjvM0KFDLUD36dMnGMxpAHX9qaeeaiV+gL0di4J5PFcUzNtxtdY3JgXz+mjfuJ5DwPy6664zw4YNM9iEP\/zwQ\/vfL774orUlrrHGGlYNPNNMM1maSW3zxx13nFlqqaWsKpJv4+zMX3zxhX1SknLWWWeZMWPGWOmX\/hZffHGDivfvf\/+7\/R3bMj4BpMp85JFHbP2nnnrK9OrVy2y33XZmq622Mj179uzEV+y9Dz\/8sLWjUrdHjx42Beeee+5pFlxwwURVeXRx\/PTTT+boo482X375ZaaNVSqZ0wd0P+igg0yc\/8IHH3xg0LAMHz7cfPbZZ5aGm2++eew8aWvUqFHm3HPPtfOFhzzTOf\/885ttttnGYFKJ2uXRHFx++eXmlltusfWnnnpqs95665mddtrJTDvttMZpIj766KMOckB7p25PAvm4dldZZRVr0kEl7Rfa+Mc\/\/mHOPvtsy1PW36uvvmrHstlmm1m+TjnllJl7Na9kHgfEbr6sPfwDTj75ZGsKYb2jNWGdLrvssh1rJgnM0bSwvp5++mkzZMgQ6xfButPSbAoomDebvy2dXSiY80gE4HvvvfdaNTIH7+OPP26eeOIJayc+\/\/zzzWyzzWZ4UOLWW28177\/\/\/jjzA4xvuukmm5\/6vPPOM4suuqg5\/PDDDVItwLHIIot0+gaABWz22Wcfe3C7A3XOOec0L7\/8spljjjnMkksuaZ3JRo4caf+bA5WxMEbABYe01157zV4aDjjggI63prFxcxm47LLLTO\/evS2gUe6++24D4OQ5aJkXIEc+bsaZVqRgzkXjsMMOMzfccIO9PK211lodzfJYx3777WdV8L\/73e8sKHO5AqihKbyabrrpOupDG3cpYozTTDONvQAwV\/7ea6+9Or3D\/frrr5tdd93V8hCg4uIEvWkfrcM555xjfS7gJZcp6mP\/n2GGGcwGG2xg\/44D82eeecbsv\/\/+tl03bvKV\/+1vfzMTTzyxvYD5+cxp45577jH9+\/e3wMnY6dfxlLqnn356JqCXCebsAejP5Yn1xiWJ8f\/666\/mtNNOM2uuuaalexyYwy\/mzwUlz\/pq6UGhnVVCAQXzSsjaPRsNBXMkRBznkKDnm28+S0QkUtTASHxI2ptuumkicf0XpfyD7L777jM77LCDBWxUyU5K4XA86aSTrKc0NmP6dGDOoYgkvPHGG5vxxhvP+MAHgAPm\/fr1s2MZPXq0leSR8n0VuJN6AUns3E7Cc5ITNmv6jV4w4iYISAHm9JtVPwvMf\/nlFztm+uYPQApgIQ1TuGhAJyRU6INUC82gF8B\/yCGHWLpwSULadtqNjz\/+2ILw3HPP3TEFR5tvv\/3Waie41MAnAOfRRx9NbB+p+NBDD7WXsiQJPPrvjJtLw0svvWQvG7\/\/\/e87eO1Afvzxx++4FDJI2uDSB2hy8XLaH8bIeuQywiUwy\/+jTDCPrj3GyUUDnnCBdN7vUTDHmXHw4MHm9ttvVyDvhkewgnk3ZHpVU5aqv+l\/oYUW6hTC44APVWjUUxspA\/VuWiibDzTrrLOOBWIX7gOg8C0qcCQzB6pO2kUSdAekA3OkToATCdMVN8Ytt9zSSrSokd2Fg\/6QIJ0K2AHcW2+91XFR8OnuNAKAli\/NJ\/EG9S+HdHRMcfUloWl8x\/jRFgDOSLquEO4FKEY1DfzuLjWAi7sAobUA+AAapHlfpevWBPN1tAFYMUkgaXMhcHSkfUAUUwbAxMtgXDCkYI72hUsVGgLG77dL24AyYMclwWk3HJjHXRQdv7MukbRdJpjHrb249n0wZ06sYXjnX0Kr2uvabvtRQMG8\/XjSZUeUJzRt5pln7iStuoMT6YhD3i8OYNdee+3EuHSn5kW6jqqAke45tK+99tpOXttcEgAV+qRtiusL1W80\/jftcI8CDippbO9IdLQTtaUjRaItQPIEoNPssk6S5dIhAf6k0LRvvvnGSpr\/\/Oc\/7dhQc\/uXFf9icsUVV1jwiz6xSR3oePDBB1u6RXnFGqD\/d955x6rl\/+\/\/\/s+qjAFlB+ZpvI5b\/FIwZzxc1pLGjare8eTYY4+1WgXaRmNw1VVXmSWWWKJT9+4S2WowR83uxucG5PYWKncXx+74TM4ALgCYJDBBoJVq57wBXfaAa\/OBK5i3OYO60vDKULPHhbRlgbmzv0KrqJrX0c8BN5Ib6sqff\/7ZHuQ4PTn1rw\/mcReHPGDugAANBICIqt4vSJ4jRoww\/\/73vzOTnkhD0lz7aWp2NAbOhyAuvtzxEDsywBAFe\/rAWQ3g8DUl\/Bsgwpy4PFGQjPF1+Prrr+3\/+45rqLal4YsSMJ911lntvLC504+v6nd04RKDtoGLJKabSSedNFHq55u6wDxu7aWBOZoRaI0Jgcuh1HTTlc4WHWs2BRTMs2mkNYQUqAPMfftz1LnJH3bUax11MWCErdSXdtMuDkXAPIt0vnd2Ul1pSJoEzB0YY39nrlEnqTymEgfmOJjR3nvvvWfNIajb55lnHjPFFFPYSxNqc1\/N7lTbdYE5zo1OW5IW9tZVwBwfAUwLCyywgDUv4AmPlgKJXUv3oYCCeffhdeUzbTWYO4c3soNFvaXjJotjHRIhkgvfYp9HivOdm8oCc6dmRw2NPbloaFCekDQpmFMPp0Dmj3ofzQTe6hTXHxJ2koQbpS32fOzb2N6xRWfZzNPU9HF8k0jmpHPNUrM7nuAY56vZk1LCdhUwR9PA5QTVOvNifR9xxBHWhFR03VV+WGgHpVNAwbx0knbfBlsJ5oAO4VT8ATCPPPJIqzZNK+4w33rrra16+9lnn+3kEMe3ZYG5c65z4VCSWOW4secJScsD5mgm8LC\/5pprrL+A7x\/gNBB4uDtfgjS6pknaSO2otunPqdmd81+cAxx8BZBQHTMmQFoK5lIHOD9BThMkcz+dq9OSsL65pKEh0dI9KKBg3j343JJZthLMXagO8cFRh7ekyTpHMpyIAJcttthinJjtssAc73ocsbhkoALFTu+ckvjtoYcesuFZ2NOpk5THPk9IWh4wp66LjcdU4ccvY\/9Gbf799993ChPkG+zfeIRjmz7xxBPNCiusYBOtIJkj6QOUzovc1SUpjG9OwBEPr3cc46IhZPwbCU8IvyNpirNrX3nllePkCSgSmsYccDhE1U5pGpgzJ3cZw+fBj+poySGgndRGAQXz2kjfvI7zeLMze9+j3R1AEgc4MpMBNnhMr7vuumb66aePJSbe34Sp4RTkigtPwrvahVb5H5cF5rTJ5QGQu\/766238\/PLLL2+92p988smOzHGkZ43zGHdjyhOSlhfM\/QsHUQCAnIuzJjqAuHwAmXHzu58EZpNNNrEhXoAt4XeOH3iE44fw\/+2dSw6DMAxEu+GoXIHTcNjqIVmy0jYUCg1Bjy0hOBOrk\/ozUIRHKx2HLdZOBCDvbRaNKcVd8ItcyBi+QZU3gj2QFIeDPaIxuW\/+DDJH+Iaiv+xz2b8Q3Jmm6UHfffmRlJrv1Qrgyg+txGEJ4RuiG+DldX8EJPP77\/HfVrileAqjcq\/5FjJ\/J\/P5bpFl+JgxtdazfP\/XavawhwgA7WAcHAjzE0aGMJmfYjGI69MVeFJIlvva1zZ0TTQmPx9CK+SHy7qDkGdFkQ8iH4ZhUegjFwsB5\/YnxkKUVMFTqQ+JMw61Mv5t04aX+7uxIWRXyaEzPxKqKLCN47iQdVwUL5JOoX2MuaP17Cg51yNz5rXPzrKe8En29iwy5z1EdFA25KAbyolrfuP9vhGQzPveP63fiACqa\/zIQSzf5IM3Tu9wERABEWiCgGTeBHZf2gIB\/iVTXEXhG1XtWfWshT2+UwREQASOQkAyPwpJ57ksAoRR53le8rjoh5tHvOxWaZgIiMBOBCTzncD5WD8IRP4Qi8kLk8stdbv7WY2WioAIiMArApK5XiECIiACIiACnSPwBGx7whSUDJcxAAAAAElFTkSuQmCC","height":300,"width":499}}
%---
