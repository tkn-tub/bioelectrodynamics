%[text] Description:
%[text] To run the Cobra toolbox you need first to download the repository as indicated [at this link](https://opencobra.github.io/cobratoolbox/stable/installation.html). Your cobra folder should be included within the file code.
clc;
clear all;
close all;

%Paremeters
fontsize=15;%for plotting

%including all parent folders up to the file 'code', this makes visible all
%files within this code, which includes the COBRA toolbox you need first to
%download as follows from the intr
addingPathParentFolderByName('code'); %[output:986f0ab3]

initCobraToolbox(false) % false, as we don't want to update %[output:6d612b06]
%This lists all the versions of functions named initCobraToolbox in your matlab installation, after the current path, the version used is the first in this list.
which -all initCobraToolbox %[output:0bb3a482]
%[text] #### **Loading the model**
% Specify the MAT file you're looking for
% matFileName = 'model_aux.mat';
% matFileName = 'iMM904.mat';
matFileName = 'model_yeastGEM_v9_0_2.mat';
%place the current Matlab folder with the project files
cd(fileparts(which(matFileName)));
%load the model
model=readCbModel(matFileName); %[output:16e16ec5]

S = model.S;
[nMets, nRxns] = size(S);
fprintf('Counter search on iMM904: %d mets, %d rxns.\n', nMets, nRxns); %[output:6f744176]
%%
%[text] #### Restrictions
%Differentiator
MAX_PRODUCERS_MD   = inf;     % m_d single-producer
MAX_PRODUCERS_MDF  = inf;     % m_Df single-producer (the new intermediate)
MAX_PRODUCERS_MDS  = inf;   % m_Ds relaxed (any producers)
EXCLUDE_CURRENCY   = true;
REQUIRE_GPR        = true;
MAX_HITS           = inf;

%Counter
MAX_PRODUCERS_MC   = inf;     % m_c controllability: single-producer in iMM904
MAX_PRODUCERS_MCS  = inf;   % m_Cs is a buffered cofactor; many producers OK
ALLOW_MD_AS_MCS    = false; % m_d ≠ m_Cs (different roles)

%NN
MAX_PRODUCERS_MNC  = inf;   % cap on producers of the convergence metabolite m_Nc (inf = no cap)
REQUIRE_GPR_NN      = true; % enzyme-mediated reactions only
EXCLUDE_CURRENCY_NN  = true; % reuse isBadMet computed in the master file
MAX_HITS_NN          = 500; % safety cap on number of motifs collected per stage
%[text] #### Currency-met filter
%Run to find names
allNames = unique(string(model.metNames));
allNames(contains(lower(allNames), ["h+","proton","hydrogen"])) %[output:2b5d9fae]
allNames(contains(lower(allNames), ["atp"])) %[output:937c7c33]

badMetNames = ["H+","H2O","ATP","ADP","amp","pi","ppi", ...
               "nad","nadhinf","nadp","nadph", ...
               "co2","o2","nh3","nh4","coa", ...
               "fad","fadh","fadh2", "oxygen"];

isBadMet = false(numel(model.mets),1);
if EXCLUDE_CURRENCY
    metNamesLower = lower(strtrim(string(model.metNames)));
    for k = 1:numel(badMetNames)
        isBadMet = isBadMet | (metNamesLower == lower(strtrim(badMetNames(k))));
    end
end


%[text] #### GPR flag
hasGPR = false(nRxns,1);
for r = 1:nRxns
    rule = strtrim(string(model.grRules{r}));
    if strlength(rule) > 0 && ~strcmpi(rule, "s0001")
        hasGPR(r) = true;
    end
end
fprintf('Enzyme-mediated rxns: %d out of %d\n', sum(hasGPR), nRxns); %[output:8934d6ad]
%[text] #### Per-metabolite producer counts
nProducers = full(sum(S > 0, 2));
%[text] #### Substrate and product sets (for shared-input check later)
substrateSet = cell(nRxns,1);
productSet   = cell(nRxns,1);
for r = 1:nRxns
    substrateSet{r} = find(S(:,r) < 0);
    productSet{r}   = find(S(:,r) > 0);
end

%[text] #### Solutions for the differentiator
Block_1_Differentiator_Mapping_anchorRDi_4rxn %[output:391715b2]

%[text] #### Solutions for the Counter
Block_2_Counter_Mapping

%[text] #### Solutions for NN
Block_3_NN_Mapping

%%
%[text] #### Printing a list of results
load Block_3_nn_full_motifs.mat
%%\
%\[text] #### Printing a list of results
%\[text] ## =========================================================
%\[text] Print the FULL combination (differentiator + counter + NN) for a
%\[text] handful of example motifs, RANKED from least to most genetically
%\[text] complex (fewest distinct genes spanned across the whole chain
%\[text] first). Every metabolite and every reaction in the chain is
%\[text] printed, with IDs, names, formulas, GPR rules, and the chemical
%\[text] equation where possible.
%\[text] =========================================================
N_PRINT = 5;                       % how many example motifs to print in full
MAX_CANDIDATES_FOR_RANKING = 2000; % safety cap on how many sourceNN entries to score before sorting

nF = numel(NN_Convergent);

fprintf('\n\n############################################################\n'); %[output:82491b5b]
fprintf('FULL DIFFERENTIATOR + COUNTER + NN MOTIFS\n'); %[output:193099ec]
fprintf('Showing %d simplest of %d scored (ranked by total distinct genes, ascending)\n', ... %[output:group:0a3879ca] %[output:3a8cd22e]
        min(N_PRINT, nF), nF); %[output:group:0a3879ca] %[output:3a8cd22e]
fprintf('############################################################\n'); %[output:2be4d4c6]

% Ordered list of (metabolite field, label) and (reaction field, label)
% pairs to print. Fields not present for a given row (e.g. m_Nc/r_No when
% hasMerge is false) are skipped automatically.
metFields = ["m_in","m_Df","m_d","m_Ds","m_c","m_Cs","m_Nc"];
metLabels = ["m_in  (differentiator input)", ...
             "m_Df  (fast intermediate)", ...
             "m_d   (differentiator output / counter input)", ...
             "m_Ds  (slow branch product)", ...
             "m_c   (counter output / NN input)", ...
             "m_Cs  (counter cofactor)", ...
             "m_Nc  (NN convergence output)"];
rxnFields = ["r_DF","r_DS","r_Dp","r_Di","r_Cp","r_Cg","r_Cr","r_Na","r_Nb","r_No"];
rxnLabels = ["r_DF (differentiator: fast)", ...
             "r_DS (differentiator: slow)", ...
             "r_Dp (differentiator: producer)", ...
             "r_Di (differentiator: integrator)", ...
             "r_Cp (counter: producer)", ...
             "r_Cg (counter: gate)", ...
             "r_Cr (counter: reset)", ...
             "r_Na (NN: neuron a)", ...
             "r_Nb (NN: neuron b)", ...
             "r_No (NN: merge/output)"];

for i = 1: nF %[output:group:30832ef2]
    rec = fullJoinRows(i);
    fprintf('\n============================================================\n'); %[output:90480593] %[output:87d38746] %[output:72e5a225] %[output:92c5e06b] %[output:1e0be3a5] %[output:8b427032] %[output:794d8672] %[output:5c949d59] %[output:6b91d1ff] %[output:981320ab] %[output:441b3eb8] %[output:762ff6e3] %[output:34d2911b] %[output:3aefed9b] %[output:20d37575] %[output:9b29afdd] %[output:50c8e3c7] %[output:3802177c] %[output:3e2821ca] %[output:69ce1586] %[output:12c73b99] %[output:43bdd8ca] %[output:7d4c6c8e] %[output:90981ff5] %[output:9a784dc8] %[output:148ad159] %[output:5fb059d9] %[output:6be7e79a] %[output:44807ecc] %[output:0c743432] %[output:7b619fd0] %[output:74fa3691] %[output:4c628799] %[output:22a8bfc7] %[output:3dac186d] %[output:8f6b2692] %[output:7fd78e3f] %[output:578ffdb6] %[output:3c13c830] %[output:722d9ba6] %[output:99f6fe94] %[output:1c3157c9] %[output:94e2264f] %[output:77b7c0f3] %[output:0e1e59fe] %[output:3b45ba51] %[output:1acc2978] %[output:5512b39e]
    fprintf(' MOTIF %d  (diff_idx=%d, cnt_idx=%d, total genes=%d)\n', ... %[output:2f5ab102] %[output:8cc72d13] %[output:56951bb7] %[output:5c73fc0a] %[output:1cabc945] %[output:7afdbf40] %[output:96d3faf8] %[output:9402ba8b] %[output:49026482] %[output:6d41f771] %[output:25b6f40b] %[output:1165b0b4] %[output:96147fd6] %[output:14b2fe30] %[output:6ec3eea5] %[output:6a7010e6] %[output:60a6c84c] %[output:4b1d3df3] %[output:465b8e3a] %[output:6c5b0824] %[output:9ef7df01] %[output:50a650c0] %[output:63cafd0e] %[output:6ecd9066] %[output:5023a33e] %[output:971bbb8d] %[output:3725d90c] %[output:5d4722e5] %[output:819564cf] %[output:7dd8c16a] %[output:57561bf6] %[output:4227a215] %[output:773f39b2] %[output:4573bb35] %[output:732834f8] %[output:91c490cb] %[output:236518b4] %[output:5d71f650] %[output:10685e5c] %[output:02e5eb98] %[output:574176bd] %[output:883afde2] %[output:1f4d246d] %[output:995a343b] %[output:657a287f] %[output:72f404a4] %[output:79224839] %[output:469aff16]
            i, rec.diff_idx, rec.cnt_idx, rec.nGenes); %[output:2f5ab102] %[output:8cc72d13] %[output:56951bb7] %[output:5c73fc0a] %[output:1cabc945] %[output:7afdbf40] %[output:96d3faf8] %[output:9402ba8b] %[output:49026482] %[output:6d41f771] %[output:25b6f40b] %[output:1165b0b4] %[output:96147fd6] %[output:14b2fe30] %[output:6ec3eea5] %[output:6a7010e6] %[output:60a6c84c] %[output:4b1d3df3] %[output:465b8e3a] %[output:6c5b0824] %[output:9ef7df01] %[output:50a650c0] %[output:63cafd0e] %[output:6ecd9066] %[output:5023a33e] %[output:971bbb8d] %[output:3725d90c] %[output:5d4722e5] %[output:819564cf] %[output:7dd8c16a] %[output:57561bf6] %[output:4227a215] %[output:773f39b2] %[output:4573bb35] %[output:732834f8] %[output:91c490cb] %[output:236518b4] %[output:5d71f650] %[output:10685e5c] %[output:02e5eb98] %[output:574176bd] %[output:883afde2] %[output:1f4d246d] %[output:995a343b] %[output:657a287f] %[output:72f404a4] %[output:79224839] %[output:469aff16]
    fprintf('============================================================\n'); %[output:31d9fd73] %[output:9733901e] %[output:45ff7400] %[output:96d34659] %[output:14b2ac8e] %[output:4bfebf99] %[output:51fba45e] %[output:346dd419] %[output:057708ae] %[output:0a8a56ae] %[output:23c8360f] %[output:6d6749b8] %[output:3897453e] %[output:1e78f4da] %[output:48868e7b] %[output:8008f20d] %[output:00f8fdba] %[output:450e25c1] %[output:5d5596df] %[output:33cdc236] %[output:1ab575f4] %[output:48b34426] %[output:51c89d90] %[output:82bc4190] %[output:9e286bcd] %[output:124e295c] %[output:18e27761] %[output:081a0053] %[output:86689416] %[output:63cb3ce2] %[output:2bfd9fc8] %[output:0682ab7b] %[output:61b974e3] %[output:45b123e5] %[output:279671dc] %[output:20f8b587] %[output:6175d1a8] %[output:17edb768] %[output:8356dff8] %[output:5af7bd9c] %[output:461543db] %[output:018c4b56] %[output:8bc7e3e3] %[output:0a514dff] %[output:830360bd] %[output:2d764a1e] %[output:94afa365] %[output:73cabbd7]

    fprintf('\n--- Metabolites ---\n'); %[output:2a57d334] %[output:8fcbaeb6] %[output:872455c6] %[output:23060af2] %[output:694e7e75] %[output:948d54b0] %[output:08bf6b4d] %[output:6c0c4b6b] %[output:4bcd9333] %[output:6cf04fc7] %[output:6cb8ddac] %[output:767e01dc] %[output:640e8294] %[output:2c44f614] %[output:150151c8] %[output:571d2d48] %[output:38bea1d6] %[output:51193ce2] %[output:62e7b090] %[output:11c2b8d8] %[output:2e868d11] %[output:48a070de] %[output:285983e7] %[output:8227a011] %[output:279d5232] %[output:47c2e24a] %[output:4a2d96df] %[output:2d76f749] %[output:7b216038] %[output:2791b88e] %[output:49d3fa21] %[output:189f0744] %[output:2aa02b52] %[output:84ebb599] %[output:430a6195] %[output:968353fa] %[output:9350e72d] %[output:102b482d] %[output:8ee2eab6] %[output:6ac165b3] %[output:79160696] %[output:4dba2507] %[output:6278097a] %[output:5ef7d6e3] %[output:8da31ad0] %[output:73a5baed] %[output:02807010] %[output:85721dcd]
    for f = 1:numel(metFields)
        mIdx = rec.(metFields(f));
        if isnan(mIdx), continue; end
        fprintf('  %-48s | %-30s | %s\n', metLabels(f), ... %[output:6fd3818a] %[output:619957f0] %[output:1cf4ac40] %[output:805f5f21] %[output:2b43278c] %[output:707dc87e] %[output:973a0b24] %[output:90aa1049] %[output:6cd288f9] %[output:60d9e700] %[output:7b39c395] %[output:33779621] %[output:6a3604ed] %[output:6e39a894] %[output:3c3ecbb7] %[output:9bb6c6ce] %[output:078bed90] %[output:50986a14] %[output:41ffd83f] %[output:429405d7] %[output:54482064] %[output:23afbc94] %[output:282295bd] %[output:862c6a68] %[output:11cbaebe] %[output:71a7d17e] %[output:5a93cbf3] %[output:1265d684] %[output:710f1a63] %[output:34854dae] %[output:24b77572] %[output:2beb4b35] %[output:807e368f] %[output:0479c72f] %[output:8b96967e] %[output:2dcb0013] %[output:3b187a08] %[output:065f1fdd] %[output:77c88e0d] %[output:3963f35b] %[output:37a54a9e] %[output:792c2b91] %[output:92595dc1] %[output:9039aeab] %[output:5840ad9d] %[output:00d557cf] %[output:5a9a9f43] %[output:873a3aa9]
                model.metNames{mIdx}, ... %[output:6fd3818a] %[output:619957f0] %[output:1cf4ac40] %[output:805f5f21] %[output:2b43278c] %[output:707dc87e] %[output:973a0b24] %[output:90aa1049] %[output:6cd288f9] %[output:60d9e700] %[output:7b39c395] %[output:33779621] %[output:6a3604ed] %[output:6e39a894] %[output:3c3ecbb7] %[output:9bb6c6ce] %[output:078bed90] %[output:50986a14] %[output:41ffd83f] %[output:429405d7] %[output:54482064] %[output:23afbc94] %[output:282295bd] %[output:862c6a68] %[output:11cbaebe] %[output:71a7d17e] %[output:5a93cbf3] %[output:1265d684] %[output:710f1a63] %[output:34854dae] %[output:24b77572] %[output:2beb4b35] %[output:807e368f] %[output:0479c72f] %[output:8b96967e] %[output:2dcb0013] %[output:3b187a08] %[output:065f1fdd] %[output:77c88e0d] %[output:3963f35b] %[output:37a54a9e] %[output:792c2b91] %[output:92595dc1] %[output:9039aeab] %[output:5840ad9d] %[output:00d557cf] %[output:5a9a9f43] %[output:873a3aa9]
                safeMetFormula(model, mIdx)); %[output:6fd3818a] %[output:619957f0] %[output:1cf4ac40] %[output:805f5f21] %[output:2b43278c] %[output:707dc87e] %[output:973a0b24] %[output:90aa1049] %[output:6cd288f9] %[output:60d9e700] %[output:7b39c395] %[output:33779621] %[output:6a3604ed] %[output:6e39a894] %[output:3c3ecbb7] %[output:9bb6c6ce] %[output:078bed90] %[output:50986a14] %[output:41ffd83f] %[output:429405d7] %[output:54482064] %[output:23afbc94] %[output:282295bd] %[output:862c6a68] %[output:11cbaebe] %[output:71a7d17e] %[output:5a93cbf3] %[output:1265d684] %[output:710f1a63] %[output:34854dae] %[output:24b77572] %[output:2beb4b35] %[output:807e368f] %[output:0479c72f] %[output:8b96967e] %[output:2dcb0013] %[output:3b187a08] %[output:065f1fdd] %[output:77c88e0d] %[output:3963f35b] %[output:37a54a9e] %[output:792c2b91] %[output:92595dc1] %[output:9039aeab] %[output:5840ad9d] %[output:00d557cf] %[output:5a9a9f43] %[output:873a3aa9]
    end

    fprintf('\n--- Reactions ---\n'); %[output:99e6fccd] %[output:03c1557f] %[output:50adb809] %[output:2214a907] %[output:0879dd9f] %[output:996c34ce] %[output:37a5c6f0] %[output:758609e1] %[output:10dd463b] %[output:17fe2835] %[output:3ead8dbf] %[output:48d3432a] %[output:856c7768] %[output:8361fbbc] %[output:57bb0fae] %[output:486272e7] %[output:7efec691] %[output:371021aa] %[output:108242b7] %[output:8a0461ad] %[output:793bfa7c] %[output:69737e47] %[output:3d7f97ac] %[output:80512382] %[output:5c95cfb2] %[output:3d22e6f9] %[output:29d4a32c] %[output:3458616f] %[output:8bf9456d] %[output:33b49a8b] %[output:5f38dcd7] %[output:511d1ab9] %[output:80a42bb2] %[output:9463d9f6] %[output:80efb81f] %[output:0003fe7a] %[output:2b59ddc7] %[output:2d907fe9] %[output:9d63ae74] %[output:5b4c6133] %[output:44ffcf62] %[output:41eac6d3] %[output:3402fc27] %[output:6f407383] %[output:8bd0324a] %[output:546a17f8] %[output:8a63e90b] %[output:65533ffa]
    for f = 1:numel(rxnFields)
        rIdx = rec.(rxnFields(f));
        if isnan(rIdx), continue; end
        fprintf('  %-32s : %-14s | %s\n', rxnLabels(f), ... %[output:964fb803] %[output:2fa6c30a] %[output:30f0434f] %[output:73591e18] %[output:8bb958de] %[output:9014db65] %[output:93d5880c] %[output:2f733ce7] %[output:86460037] %[output:965bbc39] %[output:657192d9] %[output:9c1576bd] %[output:62a42c84] %[output:79f2e308] %[output:7f48dbcb] %[output:2bbb1da1] %[output:91575452] %[output:2f9e96e5] %[output:619b6231] %[output:06319b5f] %[output:8b582948] %[output:7681edb9] %[output:4c48355e] %[output:64c521c2] %[output:5a005d6a] %[output:1ce80a59] %[output:648b807c] %[output:95583f69] %[output:0af20fff] %[output:51a8354f] %[output:6cefafbd] %[output:654dbb6a] %[output:0ea0a0a7] %[output:24a88e5f] %[output:620da475] %[output:6fa8478e] %[output:5075e836] %[output:900ccf5f] %[output:818abd41] %[output:837c2821] %[output:8dd47c55] %[output:61c25090] %[output:109ff7bc] %[output:1155f6ef] %[output:67cfd80a] %[output:7a190c11] %[output:2e5b225f] %[output:6e96a0b9] %[output:0377cd9e] %[output:8a647cb4] %[output:044b82c6] %[output:3cdc39d2] %[output:51d276c2] %[output:6caf7383] %[output:73d7674c] %[output:6b9a8477] %[output:9ef2949a] %[output:576a001a] %[output:208ceabf] %[output:31364dae] %[output:9c9de3bb] %[output:251232ef] %[output:052a4dc1] %[output:3216d79b] %[output:20840e07] %[output:356a1c5d] %[output:22066f68] %[output:660829cf] %[output:8ca8d41f] %[output:5d38b478] %[output:5c16f8ed] %[output:0f392906] %[output:0b74322e] %[output:5bbd6df2] %[output:10766e25] %[output:97509e88] %[output:7edff9fd] %[output:97e0b1a4] %[output:6be415fa] %[output:6723cec2] %[output:60638f3f] %[output:94fd2bb3] %[output:01f0a659] %[output:2a39eba0] %[output:7f7675e6] %[output:6dd39b5c] %[output:8cf84f45] %[output:78380565] %[output:1cfd6e48] %[output:54cf6d47] %[output:3ce4967e] %[output:8cf99761] %[output:652b04c4] %[output:109503f7] %[output:288558c7] %[output:2c2f7a22] %[output:124de497] %[output:005dda0f] %[output:562d04da] %[output:59132e44] %[output:1510e382] %[output:8bbc8691] %[output:2785b065] %[output:987e5db1] %[output:592173d2] %[output:6ea134fc] %[output:20bb04af] %[output:7e79642e] %[output:0db64fd5] %[output:7b2d376c] %[output:7f54a3b5] %[output:5a1ee8ff] %[output:34f8cc25] %[output:7d274839] %[output:53cd97f9] %[output:71598952] %[output:7fbba03d] %[output:1c671c66] %[output:1325f9a4] %[output:28ecd2a6] %[output:8f0fa6a8] %[output:92d1bcdb] %[output:5f01f85a] %[output:65193b44] %[output:1e09b7a7] %[output:2b2cb234] %[output:18a9961c] %[output:5a86ac05] %[output:7b789579] %[output:9536904f] %[output:6c459ae6] %[output:18723b21] %[output:14292986] %[output:97448a92] %[output:47700dae] %[output:482bbb64] %[output:3aba06ee] %[output:1264200e] %[output:8465fe82] %[output:941ddc56] %[output:8070cf7e] %[output:156b3164] %[output:308660f7] %[output:8453bb16] %[output:25950ec6] %[output:36b2fdcd] %[output:09578b1c] %[output:17f95153] %[output:756cb904] %[output:690a08c6] %[output:7cd56229] %[output:7c6daa0f] %[output:888387ac] %[output:055cfacb] %[output:141114fb] %[output:64b5f075] %[output:772e76d7] %[output:60028875] %[output:7a3d0a7c] %[output:60116d6e] %[output:4d561a93] %[output:4adef79a] %[output:2c74ca67] %[output:825896b9] %[output:2bfe58f1] %[output:96f30d5d] %[output:28327f8d] %[output:5b17f74e] %[output:9e82fe58] %[output:7229a42b] %[output:0e9e6037] %[output:8bc3b2f3] %[output:2aa4393d] %[output:8e408e64] %[output:529926c4] %[output:56ef8566] %[output:510adde2] %[output:2619dac1] %[output:59de3bd1] %[output:20e19cdc] %[output:491bd3f1] %[output:877bbda3] %[output:83bea097] %[output:225fa24c] %[output:6a8f15d3] %[output:723fcdbf] %[output:29517250] %[output:415f3196] %[output:7829b937] %[output:4dae2b85] %[output:8f5fd96c] %[output:811af5b5] %[output:042758ba] %[output:0a1727a4] %[output:3c90121d] %[output:9431c610] %[output:18ed485d] %[output:1b7466a0] %[output:045f3b73] %[output:5dd6406b] %[output:09661b1f] %[output:63346b63] %[output:2cf2746d] %[output:0ad804bd] %[output:68c00048] %[output:547544ba] %[output:7337ec59] %[output:2c713c3d] %[output:1f790cba] %[output:7ff0ca6e] %[output:59d62410] %[output:624dce73] %[output:100a3259] %[output:6cc20763] %[output:1e475c9d] %[output:7f69e7a8] %[output:0ac58742] %[output:3239fec5] %[output:90345e89] %[output:1c2ff48c] %[output:81ca3f42] %[output:7c55fbae] %[output:56442f23] %[output:49ed3f3a] %[output:1e6deb54] %[output:17203297] %[output:1a467ed6] %[output:438d7e38] %[output:563c053e] %[output:6776a85c] %[output:81154071] %[output:2309ebd3] %[output:0d57fddd] %[output:733cf8df] %[output:1ee35c4a] %[output:709e25ac] %[output:4dbb3e9d] %[output:7c200491] %[output:4d3a21a3] %[output:64529e4f] %[output:00aacdba] %[output:8ff21982] %[output:2243b79f] %[output:805fd59c] %[output:54723426] %[output:938f5a07] %[output:03afeb76] %[output:87fc5a05] %[output:2ea93c13] %[output:01e41403] %[output:0972ae65] %[output:787e4ca6] %[output:6cfb7f07] %[output:158b2327] %[output:82af6900] %[output:1ec93e3f] %[output:0c62363d] %[output:629ab37e] %[output:28417d3e] %[output:5eac50b9] %[output:6c12301f] %[output:4a2c42e0] %[output:54054711] %[output:70c99a1e] %[output:82c0d72f] %[output:23fe02c5] %[output:13ca59f2] %[output:9b88ed30] %[output:74955944] %[output:42ad7701] %[output:39591b81] %[output:69cdad03] %[output:9b897334] %[output:3140324a] %[output:3bbe9e20] %[output:5f665678] %[output:4c483888] %[output:462430cc] %[output:8039e85a] %[output:342fdf63] %[output:32e401f6] %[output:7c109e76] %[output:02118ba6] %[output:835d07a6] %[output:783dfb41] %[output:1ee3b958] %[output:8bf81094] %[output:5a303777] %[output:4b75ddbe] %[output:73055ee5] %[output:19948de7] %[output:75beb517] %[output:70b7bcc0] %[output:31f81824] %[output:9cebd55f] %[output:957104f1] %[output:14305eaa] %[output:61054bab] %[output:0475ad07] %[output:649dc9e9] %[output:5e5e25b1] %[output:958144ad] %[output:2b5887ca] %[output:1eaa5ed2] %[output:17fe234c] %[output:48c595be] %[output:8382b765] %[output:9b22f22c] %[output:48799ffa] %[output:892431f2] %[output:0c332368] %[output:60f0a685] %[output:184d77d5] %[output:35fa9a7b] %[output:9ff3da68] %[output:456c5c69] %[output:34007a44] %[output:58cccef7] %[output:6357edde] %[output:15acaac5] %[output:70706130] %[output:8abd0b2c] %[output:0ed45ef0] %[output:00f07ee0] %[output:52b7a64d] %[output:17f40f08] %[output:1623b9d5] %[output:91ce738c] %[output:174c1948] %[output:6922ac8b] %[output:0baccec1] %[output:8550f37b] %[output:2aeda85d] %[output:01058bcc] %[output:4fb4b327] %[output:19903fdb] %[output:2894550b] %[output:8514f1bb] %[output:6721459e] %[output:32663a74] %[output:5542d5da] %[output:614d4d7f] %[output:7f773905] %[output:3cd9242c] %[output:3eb168dd] %[output:6ac5d2fd] %[output:20afdda3] %[output:2d28c47a] %[output:3e9a918a] %[output:44ba9a0e] %[output:25cab39a] %[output:96f00851] %[output:8c5ab970] %[output:5bd591aa] %[output:8f89d508] %[output:1a285fde] %[output:1976301d] %[output:9a8397da] %[output:393c1bb0] %[output:24fbfaba] %[output:7968c866] %[output:37f867fa] %[output:2ae8b8e9] %[output:562c2813] %[output:0a97bc81] %[output:27701269] %[output:06bfb699] %[output:484a2858] %[output:2d67fd90] %[output:8c7ff89d] %[output:8ba440b2] %[output:87e3c15b] %[output:8d7afd63] %[output:5dd21c1e] %[output:9f1d5ebe] %[output:3e07d1a3] %[output:95e20770] %[output:1c52fc70] %[output:56b44529] %[output:18b1c781] %[output:219952e6] %[output:248be4f5] %[output:64b95e3d] %[output:7a3a74e7] %[output:7c8fb0c0] %[output:4a1d01d9] %[output:31ededcf] %[output:99a1bddd] %[output:86d61545] %[output:84603a2c] %[output:10278aa8] %[output:41aedba7] %[output:6d075e35] %[output:11b755ab] %[output:56175571] %[output:67b68563] %[output:49a848b3] %[output:1ea8f1f9] %[output:83a146cd] %[output:1e149a03] %[output:9f929bf7] %[output:500d110d] %[output:6e2e6f98] %[output:6f44d599] %[output:1efffd2a] %[output:13b16224] %[output:235b6b80] %[output:76e886ba] %[output:2874229c] %[output:281445fb] %[output:07818ad6] %[output:4b8db7db] %[output:77db4229] %[output:2c4f9426] %[output:13535b23] %[output:1092e82c] %[output:6546d9dc] %[output:936e8057] %[output:1c9db429] %[output:280dc209] %[output:37a41491] %[output:450b4592] %[output:6cbc338b] %[output:2ff5844b] %[output:758334d7] %[output:0aaf3c01] %[output:714c5380] %[output:98c39a44] %[output:27647209] %[output:7b6eb432] %[output:16381bdf] %[output:6dc67b30] %[output:332bb081] %[output:63a01f3d] %[output:0696a6a5] %[output:96de334a] %[output:501d9f5f] %[output:4db06b73] %[output:57009c43] %[output:20581251] %[output:6a07e951] %[output:905f6834] %[output:9202b44f] %[output:2d7f544d] %[output:9f59faa4] %[output:785dd518] %[output:74945751] %[output:094181d6] %[output:220fc470] %[output:6706c0de] %[output:2dc8503b] %[output:76f6306a] %[output:3f057189] %[output:82da4b83] %[output:211e9d7a] %[output:25ab0747] %[output:7b1acf87] %[output:21452143] %[output:30fe8ec8] %[output:9af411f2] %[output:2cecce94] %[output:0878a074] %[output:7566d7a7] %[output:2d99881a] %[output:2b78ba07] %[output:02257bf8] %[output:04597955] %[output:6c5c0671] %[output:268bbce2] %[output:1c2e2eb3] %[output:3834a53c] %[output:95e395df] %[output:1b349ca6] %[output:9e3489df] %[output:9444c103] %[output:53086a1e] %[output:169dd99a] %[output:1a677094] %[output:4e25ffd2] %[output:4f0e8a9e]
                model.rxns{rIdx}, model.rxnNames{rIdx}); %[output:964fb803] %[output:2fa6c30a] %[output:30f0434f] %[output:73591e18] %[output:8bb958de] %[output:9014db65] %[output:93d5880c] %[output:2f733ce7] %[output:86460037] %[output:965bbc39] %[output:657192d9] %[output:9c1576bd] %[output:62a42c84] %[output:79f2e308] %[output:7f48dbcb] %[output:2bbb1da1] %[output:91575452] %[output:2f9e96e5] %[output:619b6231] %[output:06319b5f] %[output:8b582948] %[output:7681edb9] %[output:4c48355e] %[output:64c521c2] %[output:5a005d6a] %[output:1ce80a59] %[output:648b807c] %[output:95583f69] %[output:0af20fff] %[output:51a8354f] %[output:6cefafbd] %[output:654dbb6a] %[output:0ea0a0a7] %[output:24a88e5f] %[output:620da475] %[output:6fa8478e] %[output:5075e836] %[output:900ccf5f] %[output:818abd41] %[output:837c2821] %[output:8dd47c55] %[output:61c25090] %[output:109ff7bc] %[output:1155f6ef] %[output:67cfd80a] %[output:7a190c11] %[output:2e5b225f] %[output:6e96a0b9] %[output:0377cd9e] %[output:8a647cb4] %[output:044b82c6] %[output:3cdc39d2] %[output:51d276c2] %[output:6caf7383] %[output:73d7674c] %[output:6b9a8477] %[output:9ef2949a] %[output:576a001a] %[output:208ceabf] %[output:31364dae] %[output:9c9de3bb] %[output:251232ef] %[output:052a4dc1] %[output:3216d79b] %[output:20840e07] %[output:356a1c5d] %[output:22066f68] %[output:660829cf] %[output:8ca8d41f] %[output:5d38b478] %[output:5c16f8ed] %[output:0f392906] %[output:0b74322e] %[output:5bbd6df2] %[output:10766e25] %[output:97509e88] %[output:7edff9fd] %[output:97e0b1a4] %[output:6be415fa] %[output:6723cec2] %[output:60638f3f] %[output:94fd2bb3] %[output:01f0a659] %[output:2a39eba0] %[output:7f7675e6] %[output:6dd39b5c] %[output:8cf84f45] %[output:78380565] %[output:1cfd6e48] %[output:54cf6d47] %[output:3ce4967e] %[output:8cf99761] %[output:652b04c4] %[output:109503f7] %[output:288558c7] %[output:2c2f7a22] %[output:124de497] %[output:005dda0f] %[output:562d04da] %[output:59132e44] %[output:1510e382] %[output:8bbc8691] %[output:2785b065] %[output:987e5db1] %[output:592173d2] %[output:6ea134fc] %[output:20bb04af] %[output:7e79642e] %[output:0db64fd5] %[output:7b2d376c] %[output:7f54a3b5] %[output:5a1ee8ff] %[output:34f8cc25] %[output:7d274839] %[output:53cd97f9] %[output:71598952] %[output:7fbba03d] %[output:1c671c66] %[output:1325f9a4] %[output:28ecd2a6] %[output:8f0fa6a8] %[output:92d1bcdb] %[output:5f01f85a] %[output:65193b44] %[output:1e09b7a7] %[output:2b2cb234] %[output:18a9961c] %[output:5a86ac05] %[output:7b789579] %[output:9536904f] %[output:6c459ae6] %[output:18723b21] %[output:14292986] %[output:97448a92] %[output:47700dae] %[output:482bbb64] %[output:3aba06ee] %[output:1264200e] %[output:8465fe82] %[output:941ddc56] %[output:8070cf7e] %[output:156b3164] %[output:308660f7] %[output:8453bb16] %[output:25950ec6] %[output:36b2fdcd] %[output:09578b1c] %[output:17f95153] %[output:756cb904] %[output:690a08c6] %[output:7cd56229] %[output:7c6daa0f] %[output:888387ac] %[output:055cfacb] %[output:141114fb] %[output:64b5f075] %[output:772e76d7] %[output:60028875] %[output:7a3d0a7c] %[output:60116d6e] %[output:4d561a93] %[output:4adef79a] %[output:2c74ca67] %[output:825896b9] %[output:2bfe58f1] %[output:96f30d5d] %[output:28327f8d] %[output:5b17f74e] %[output:9e82fe58] %[output:7229a42b] %[output:0e9e6037] %[output:8bc3b2f3] %[output:2aa4393d] %[output:8e408e64] %[output:529926c4] %[output:56ef8566] %[output:510adde2] %[output:2619dac1] %[output:59de3bd1] %[output:20e19cdc] %[output:491bd3f1] %[output:877bbda3] %[output:83bea097] %[output:225fa24c] %[output:6a8f15d3] %[output:723fcdbf] %[output:29517250] %[output:415f3196] %[output:7829b937] %[output:4dae2b85] %[output:8f5fd96c] %[output:811af5b5] %[output:042758ba] %[output:0a1727a4] %[output:3c90121d] %[output:9431c610] %[output:18ed485d] %[output:1b7466a0] %[output:045f3b73] %[output:5dd6406b] %[output:09661b1f] %[output:63346b63] %[output:2cf2746d] %[output:0ad804bd] %[output:68c00048] %[output:547544ba] %[output:7337ec59] %[output:2c713c3d] %[output:1f790cba] %[output:7ff0ca6e] %[output:59d62410] %[output:624dce73] %[output:100a3259] %[output:6cc20763] %[output:1e475c9d] %[output:7f69e7a8] %[output:0ac58742] %[output:3239fec5] %[output:90345e89] %[output:1c2ff48c] %[output:81ca3f42] %[output:7c55fbae] %[output:56442f23] %[output:49ed3f3a] %[output:1e6deb54] %[output:17203297] %[output:1a467ed6] %[output:438d7e38] %[output:563c053e] %[output:6776a85c] %[output:81154071] %[output:2309ebd3] %[output:0d57fddd] %[output:733cf8df] %[output:1ee35c4a] %[output:709e25ac] %[output:4dbb3e9d] %[output:7c200491] %[output:4d3a21a3] %[output:64529e4f] %[output:00aacdba] %[output:8ff21982] %[output:2243b79f] %[output:805fd59c] %[output:54723426] %[output:938f5a07] %[output:03afeb76] %[output:87fc5a05] %[output:2ea93c13] %[output:01e41403] %[output:0972ae65] %[output:787e4ca6] %[output:6cfb7f07] %[output:158b2327] %[output:82af6900] %[output:1ec93e3f] %[output:0c62363d] %[output:629ab37e] %[output:28417d3e] %[output:5eac50b9] %[output:6c12301f] %[output:4a2c42e0] %[output:54054711] %[output:70c99a1e] %[output:82c0d72f] %[output:23fe02c5] %[output:13ca59f2] %[output:9b88ed30] %[output:74955944] %[output:42ad7701] %[output:39591b81] %[output:69cdad03] %[output:9b897334] %[output:3140324a] %[output:3bbe9e20] %[output:5f665678] %[output:4c483888] %[output:462430cc] %[output:8039e85a] %[output:342fdf63] %[output:32e401f6] %[output:7c109e76] %[output:02118ba6] %[output:835d07a6] %[output:783dfb41] %[output:1ee3b958] %[output:8bf81094] %[output:5a303777] %[output:4b75ddbe] %[output:73055ee5] %[output:19948de7] %[output:75beb517] %[output:70b7bcc0] %[output:31f81824] %[output:9cebd55f] %[output:957104f1] %[output:14305eaa] %[output:61054bab] %[output:0475ad07] %[output:649dc9e9] %[output:5e5e25b1] %[output:958144ad] %[output:2b5887ca] %[output:1eaa5ed2] %[output:17fe234c] %[output:48c595be] %[output:8382b765] %[output:9b22f22c] %[output:48799ffa] %[output:892431f2] %[output:0c332368] %[output:60f0a685] %[output:184d77d5] %[output:35fa9a7b] %[output:9ff3da68] %[output:456c5c69] %[output:34007a44] %[output:58cccef7] %[output:6357edde] %[output:15acaac5] %[output:70706130] %[output:8abd0b2c] %[output:0ed45ef0] %[output:00f07ee0] %[output:52b7a64d] %[output:17f40f08] %[output:1623b9d5] %[output:91ce738c] %[output:174c1948] %[output:6922ac8b] %[output:0baccec1] %[output:8550f37b] %[output:2aeda85d] %[output:01058bcc] %[output:4fb4b327] %[output:19903fdb] %[output:2894550b] %[output:8514f1bb] %[output:6721459e] %[output:32663a74] %[output:5542d5da] %[output:614d4d7f] %[output:7f773905] %[output:3cd9242c] %[output:3eb168dd] %[output:6ac5d2fd] %[output:20afdda3] %[output:2d28c47a] %[output:3e9a918a] %[output:44ba9a0e] %[output:25cab39a] %[output:96f00851] %[output:8c5ab970] %[output:5bd591aa] %[output:8f89d508] %[output:1a285fde] %[output:1976301d] %[output:9a8397da] %[output:393c1bb0] %[output:24fbfaba] %[output:7968c866] %[output:37f867fa] %[output:2ae8b8e9] %[output:562c2813] %[output:0a97bc81] %[output:27701269] %[output:06bfb699] %[output:484a2858] %[output:2d67fd90] %[output:8c7ff89d] %[output:8ba440b2] %[output:87e3c15b] %[output:8d7afd63] %[output:5dd21c1e] %[output:9f1d5ebe] %[output:3e07d1a3] %[output:95e20770] %[output:1c52fc70] %[output:56b44529] %[output:18b1c781] %[output:219952e6] %[output:248be4f5] %[output:64b95e3d] %[output:7a3a74e7] %[output:7c8fb0c0] %[output:4a1d01d9] %[output:31ededcf] %[output:99a1bddd] %[output:86d61545] %[output:84603a2c] %[output:10278aa8] %[output:41aedba7] %[output:6d075e35] %[output:11b755ab] %[output:56175571] %[output:67b68563] %[output:49a848b3] %[output:1ea8f1f9] %[output:83a146cd] %[output:1e149a03] %[output:9f929bf7] %[output:500d110d] %[output:6e2e6f98] %[output:6f44d599] %[output:1efffd2a] %[output:13b16224] %[output:235b6b80] %[output:76e886ba] %[output:2874229c] %[output:281445fb] %[output:07818ad6] %[output:4b8db7db] %[output:77db4229] %[output:2c4f9426] %[output:13535b23] %[output:1092e82c] %[output:6546d9dc] %[output:936e8057] %[output:1c9db429] %[output:280dc209] %[output:37a41491] %[output:450b4592] %[output:6cbc338b] %[output:2ff5844b] %[output:758334d7] %[output:0aaf3c01] %[output:714c5380] %[output:98c39a44] %[output:27647209] %[output:7b6eb432] %[output:16381bdf] %[output:6dc67b30] %[output:332bb081] %[output:63a01f3d] %[output:0696a6a5] %[output:96de334a] %[output:501d9f5f] %[output:4db06b73] %[output:57009c43] %[output:20581251] %[output:6a07e951] %[output:905f6834] %[output:9202b44f] %[output:2d7f544d] %[output:9f59faa4] %[output:785dd518] %[output:74945751] %[output:094181d6] %[output:220fc470] %[output:6706c0de] %[output:2dc8503b] %[output:76f6306a] %[output:3f057189] %[output:82da4b83] %[output:211e9d7a] %[output:25ab0747] %[output:7b1acf87] %[output:21452143] %[output:30fe8ec8] %[output:9af411f2] %[output:2cecce94] %[output:0878a074] %[output:7566d7a7] %[output:2d99881a] %[output:2b78ba07] %[output:02257bf8] %[output:04597955] %[output:6c5c0671] %[output:268bbce2] %[output:1c2e2eb3] %[output:3834a53c] %[output:95e395df] %[output:1b349ca6] %[output:9e3489df] %[output:9444c103] %[output:53086a1e] %[output:169dd99a] %[output:1a677094] %[output:4e25ffd2] %[output:4f0e8a9e]
        fprintf('      GPR: %s\n', safeGPR(model, rIdx)); %[output:579df49a] %[output:6e54ee94] %[output:59b39ee9] %[output:925d121e] %[output:6403bf32] %[output:1ae71e83] %[output:5a3d340d] %[output:4dc2e804] %[output:66364b24] %[output:6b0a14bd] %[output:2d9146c8] %[output:808bd563] %[output:22a201bf] %[output:849175bc] %[output:47f98be8] %[output:9fa608e5] %[output:70bbf2e1] %[output:2af384bc] %[output:92b262fc] %[output:687bcc2b] %[output:9dbc1e29] %[output:8b37fa4e] %[output:5a9b1eb4] %[output:852cd886] %[output:7ca31de2] %[output:9e5ebd64] %[output:972a5e47] %[output:8db61185] %[output:43084eaf] %[output:4b002ffc] %[output:3e6b709b] %[output:1dad3065] %[output:8df2beb6] %[output:9ae86943] %[output:97dbf21c] %[output:46b4e7e3] %[output:495b55aa] %[output:9786816e] %[output:92aa1e38] %[output:03440e61] %[output:64f56ee7] %[output:4d8a3c9b] %[output:9972d2fa] %[output:4125fb39] %[output:87594fe4] %[output:38c0842e] %[output:98f5e7f9] %[output:0961979e] %[output:806f2c93] %[output:96779ef6] %[output:23071a4f] %[output:00dac4ee] %[output:19163115] %[output:892cba43] %[output:6bc6e900] %[output:3067d587] %[output:59987011] %[output:7c442532] %[output:0f32864a] %[output:180aee39] %[output:04c7a9a7] %[output:678bb670] %[output:6b48ece5] %[output:68ff350f] %[output:514aab9d] %[output:9c342979] %[output:64f1fd7c] %[output:9d16ed42] %[output:38527866] %[output:04bacea8] %[output:049a7b31] %[output:20bd383f] %[output:814cfcfe] %[output:11a3ea93] %[output:32b8b06e] %[output:6e621266] %[output:77cd6f4c] %[output:3200e3f0] %[output:7b5ec38f] %[output:2be1a440] %[output:442b20b2] %[output:36e7dc5c] %[output:916ea3b1] %[output:7004f370] %[output:7c55697a] %[output:9997b74e] %[output:59d108ac] %[output:0307d386] %[output:850b854d] %[output:67010e23] %[output:4906a666] %[output:87bf322e] %[output:65af014d] %[output:5c44038e] %[output:6daf15a0] %[output:4bcc4f48] %[output:1421e4c1] %[output:1b71d38a] %[output:21d1fe90] %[output:1cfbb3ea] %[output:5a162b5c] %[output:542eb6f1] %[output:9557989e] %[output:14bd5abd] %[output:65a2ad73] %[output:7b77b78e] %[output:5cfe118a] %[output:3212ffb6] %[output:6f44788f] %[output:5a93c114] %[output:7ad2409a] %[output:98161b69] %[output:73a953e1] %[output:2a5edc41] %[output:83cabba2] %[output:6b7db9ee] %[output:299de291] %[output:627299b7] %[output:74660c87] %[output:4a71bad3] %[output:134c19ad] %[output:84ca39c9] %[output:08e869d2] %[output:0ca3fbb0] %[output:13eb50e5] %[output:63034c2d] %[output:35e6ed0d] %[output:38e22857] %[output:76d1ff4c] %[output:0ff528b6] %[output:2c6692ac] %[output:809bd8a1] %[output:2b2cc5ee] %[output:314f7d37] %[output:3f27e2be] %[output:33e12251] %[output:49c92ec8] %[output:7bcd1b99] %[output:895feaf5] %[output:2539792e] %[output:00cd6647] %[output:197774a5] %[output:6ee502fa] %[output:208e3c33] %[output:2cf198fe] %[output:3168c4a7] %[output:60eedb4d] %[output:9967f689] %[output:94b5812c] %[output:923ba3d0] %[output:8574ef95] %[output:52b6686d] %[output:247a78a0] %[output:4fda50c6] %[output:56a7931e] %[output:6dd38e0e] %[output:4c66c717] %[output:8a826441] %[output:7f5812c7] %[output:9d2222dd] %[output:77a02292] %[output:370fc088] %[output:1effe680] %[output:404dc7ff] %[output:57374d5b] %[output:2610d67b] %[output:37cec4d3] %[output:3e762204] %[output:39d69461] %[output:1840d493] %[output:997db68b] %[output:71dbe905] %[output:4208a70d] %[output:3c66d2dd] %[output:637f710d] %[output:1f9ceee8] %[output:083d58d6] %[output:6032c1a5] %[output:6590dd0c] %[output:0c56f190] %[output:5d40f319] %[output:3c271d5b] %[output:531dfbc1] %[output:6bbfb28f] %[output:9f12b557] %[output:17f12838] %[output:4a6017fe] %[output:7223afa0] %[output:2ba89321] %[output:38662d67] %[output:25f02230] %[output:17354a54] %[output:98236aac] %[output:2afff1a9] %[output:65a42c56] %[output:55200965] %[output:03cfe47a] %[output:721cc78e] %[output:422ab7f2] %[output:8487a50e] %[output:2d11f1cf] %[output:4d949db6] %[output:987a3b13] %[output:3d762783] %[output:0de969d0] %[output:27bb10c3] %[output:469a23df] %[output:0fd2714b] %[output:72920c73] %[output:279ee2a2] %[output:07119c52] %[output:8f5711ec] %[output:8f383d3d] %[output:94d778d0] %[output:9d1dea25] %[output:65f1cf8a] %[output:69014988] %[output:53423dd2] %[output:4f2224c7] %[output:5baecd8d] %[output:5b29282f] %[output:2587dbeb] %[output:95469740] %[output:601aad3e] %[output:8cc0e72f] %[output:0cb49b29] %[output:0a331101] %[output:1424d20b] %[output:0b676ae2] %[output:01369c25] %[output:67897898] %[output:0dc52fb3] %[output:5e905780] %[output:88c74c05] %[output:1378f98b] %[output:5b4fdcb8] %[output:8bf8af8e] %[output:2e00c8d8] %[output:0dfadf2c] %[output:1d939c0b] %[output:085aa13e] %[output:96e3cfde] %[output:2447c81b] %[output:8f39e3ca] %[output:2bf83424] %[output:215390cf] %[output:90c1500a] %[output:787f59b8] %[output:79f2b9ee] %[output:88379c21] %[output:7eaf6d82] %[output:8bde2a3a] %[output:3302ef6f] %[output:9774d8ff] %[output:909e7f97] %[output:51bcb6bd] %[output:02652d33] %[output:16710877] %[output:584b9ff7] %[output:97040cbf] %[output:71fb6e2a] %[output:35032c67] %[output:4060240e] %[output:6d29648c] %[output:82b4c3a6] %[output:61ebbaf8] %[output:37533561] %[output:0cba1e27] %[output:4634113f] %[output:80172205] %[output:993cb599] %[output:4a2f6fb7] %[output:3c44c574] %[output:6498b959] %[output:1a1a1bba] %[output:11967c4a] %[output:5e2f29d3] %[output:8a02a723] %[output:3fa6693a] %[output:5d7ec33c] %[output:9f52ed87] %[output:0f48f0c7] %[output:5490f5c1] %[output:1d637d4c] %[output:38160bd8] %[output:9dd636a9] %[output:5825699b] %[output:286055e9] %[output:69c52bbe] %[output:16dcb73e] %[output:84c85ab3] %[output:03cb4bab] %[output:9148aff4] %[output:3eb600a9] %[output:6ac931ad] %[output:16e93374] %[output:2a1d0815] %[output:9bfafbc5] %[output:04fa4a53] %[output:00e4d0fa] %[output:81fa684b] %[output:21c76914] %[output:69b2020d] %[output:2226baf7] %[output:083fd791] %[output:3b5caaa8] %[output:96f09537] %[output:54dcc1a3] %[output:9d757aa5] %[output:6283f531] %[output:53f01bed] %[output:47462cd9] %[output:0038c101] %[output:0cf01430] %[output:2631ccf2] %[output:597c42aa] %[output:12888b9a] %[output:61c7c548] %[output:2454b840] %[output:0bcf8625] %[output:9526e22a] %[output:50af5348] %[output:55dc1250] %[output:947d5dc3] %[output:7dcf83a7] %[output:88b2f203] %[output:806b348d] %[output:1bcd02d3] %[output:9ac5bdef] %[output:44c3cd9f] %[output:5e255b2d] %[output:0e723180] %[output:5f29262c] %[output:09af1cbc] %[output:03b0f43c] %[output:7bae529e] %[output:3fc10494] %[output:692c8754] %[output:3c6d5b6e] %[output:7558509c] %[output:2960e4eb] %[output:6f8480fc] %[output:0fe96742] %[output:73a3cc9a] %[output:3d47bc10] %[output:75ff1a96] %[output:7c9826f3] %[output:08b131b9] %[output:841da2ba] %[output:6bd51730] %[output:86662e50] %[output:906e0875] %[output:6100f547] %[output:32825530] %[output:073d68f9] %[output:8da1fb49] %[output:7b1e2d9e] %[output:93645534] %[output:69f7d401] %[output:8e9de97d] %[output:44f7c261] %[output:5404073d] %[output:3b70c5ef] %[output:0cd6c570] %[output:5e0afd21] %[output:109ba664] %[output:33a516c6] %[output:40d6cef3] %[output:24f00a60] %[output:30ae46de] %[output:172727bb] %[output:3a95dadc] %[output:63c778f6] %[output:0af3847f] %[output:6a26ccd9] %[output:6466dfef] %[output:2b09d326] %[output:0bc90a0a] %[output:05c32c83] %[output:2c6f6cb3] %[output:51648625] %[output:698e57eb] %[output:010e08a9] %[output:0ad12e5a] %[output:6d2cf2f7] %[output:1ddc4e2b] %[output:18652103] %[output:4ee824fc] %[output:88df9472] %[output:4390b76d] %[output:7618cd5e] %[output:1443a754] %[output:4d4b6fa7] %[output:9ad3b3c5] %[output:332678dc] %[output:31f28b69] %[output:020e17cf] %[output:67c95a80] %[output:277744ba] %[output:217a8be8] %[output:672fd99e] %[output:4c804c69] %[output:157565a4] %[output:5818c3f1] %[output:0ebc8bf8] %[output:1f706c8a] %[output:5ed1f4b6] %[output:2a8c1dc0] %[output:32170d8d] %[output:50458c07] %[output:6f7ebc90] %[output:16df5828] %[output:2715ee4a] %[output:265ee563] %[output:60e3950f] %[output:107ac611] %[output:3bdb9577] %[output:249bfa72] %[output:8dc9886b] %[output:09858bf1] %[output:1651321d] %[output:98258f39] %[output:188a2fda] %[output:01d6e0c0] %[output:63cf6405] %[output:33de7d25] %[output:7963f6ac] %[output:33a2e385] %[output:406096cb] %[output:31dfa6a5] %[output:7258a7e4] %[output:9347860b] %[output:305ffae5] %[output:78914f6b] %[output:77f7a612] %[output:068cedf5] %[output:79dec1f6] %[output:6c9085a5] %[output:397abd96] %[output:6e2c958f] %[output:60b994b2] %[output:01526602] %[output:9fc7c8a4] %[output:64452154] %[output:6bd6dff5] %[output:3f4cb34a] %[output:002295a1] %[output:54ebc1ff] %[output:4b643aee] %[output:5b0c963b] %[output:5ccb59cd] %[output:32acb1a7] %[output:221e87da] %[output:81e14dc1] %[output:9899458c] %[output:088e37b8] %[output:989075f4] %[output:863f713d] %[output:1020c36c] %[output:46072f23] %[output:4a72b0f8] %[output:86030b55] %[output:6a414d77] %[output:1f634206] %[output:741f71aa] %[output:4c9d1731] %[output:7032cdf1] %[output:76bd3147] %[output:747c79cb] %[output:7aea2c68] %[output:65f4b7c4] %[output:624bd7bd] %[output:405485c3] %[output:54b14d93] %[output:5daf6930] %[output:6ed27562] %[output:27859618] %[output:1fb0cd1e] %[output:1734c8e7] %[output:3e237f4a]
        fprintf('      EC : %s\n', safeECCode(model, rIdx)); %[output:04de2ea6] %[output:962adbdf] %[output:65a81a21] %[output:4ba27199] %[output:2ad6032f] %[output:228a90b3] %[output:5d753d9b] %[output:6700baa7] %[output:30d2eb08] %[output:2029ab25] %[output:2f7409df] %[output:5408d0af] %[output:99b2f8a2] %[output:2645637c] %[output:46002d9c] %[output:9209d04c] %[output:6d83ae88] %[output:93ab954d] %[output:5007c0d9] %[output:4be7c015] %[output:0de14ea5] %[output:51f71ce4] %[output:8f1ada4a] %[output:3749cea5] %[output:24e4049a] %[output:0ee78ee8] %[output:2a881c51] %[output:75f2adc8] %[output:3fac79ef] %[output:04458a67] %[output:6b36469d] %[output:08f05f0e] %[output:2b19d544] %[output:08a1a377] %[output:349b76d6] %[output:6469a447] %[output:681ba772] %[output:13590dff] %[output:0dd8ed02] %[output:9c84642a] %[output:6450ea47] %[output:388b3e5d] %[output:7d35be0a] %[output:6c4422c6] %[output:1ae29f36] %[output:3756a6f6] %[output:0c7bb582] %[output:8bc95633] %[output:2f9043e6] %[output:13ad370a] %[output:58b1fce7] %[output:6f362ef1] %[output:9dd8af12] %[output:011b3745] %[output:0a0c449f] %[output:64f71a5f] %[output:987b6de9] %[output:33d0cf16] %[output:4994efd2] %[output:5b84f492] %[output:88d7adb4] %[output:047ba959] %[output:24f7dd0f] %[output:9bb19f9b] %[output:81724c1a] %[output:6cdbdd1e] %[output:9c4e81d9] %[output:2ca73f8e] %[output:37a33a25] %[output:5a4088d1] %[output:80c159c0] %[output:0d79c1ff] %[output:7cbf001a] %[output:316a5dbe] %[output:30e40743] %[output:9557eb60] %[output:3b81fabb] %[output:50990ed8] %[output:47842898] %[output:9f173541] %[output:28292ae4] %[output:391354f6] %[output:03bd4efe] %[output:9d44bcb0] %[output:11ccd763] %[output:36265e10] %[output:3b72d254] %[output:71f17d2d] %[output:5210269f] %[output:212e233a] %[output:225af4b3] %[output:80381dff] %[output:9aefd86f] %[output:7bcfb366] %[output:55c447bc] %[output:7cab1ecf] %[output:54f1bb70] %[output:5fb97a03] %[output:8b04fa95] %[output:09c5fa8b] %[output:79004910] %[output:660c30ee] %[output:0b8cc50e] %[output:30c6d3b5] %[output:3fae533b] %[output:5c3258a3] %[output:99d6e468] %[output:16012c06] %[output:696bda48] %[output:7a3cc764] %[output:64e11a91] %[output:3e4638ac] %[output:422c3da0] %[output:6018e26b] %[output:00a9ceb4] %[output:5eb400bf] %[output:7426e05e] %[output:0597eebb] %[output:3b28b432] %[output:7beb33cc] %[output:63e4edbd] %[output:8cbcd850] %[output:39ee2afe] %[output:698fe15b] %[output:3bf4a395] %[output:76cefc60] %[output:0b61c299] %[output:0a599559] %[output:3c27f940] %[output:4d7dd10d] %[output:588e3100] %[output:9c04d3c6] %[output:233ad83c] %[output:7c100058] %[output:55e04666] %[output:91d3d3b7] %[output:324b21e9] %[output:34b4676c] %[output:5b44c6cb] %[output:92701816] %[output:9911cb42] %[output:6aa21168] %[output:0f2eed59] %[output:1f44d949] %[output:80ced352] %[output:1becb455] %[output:9218fb29] %[output:55b6f093] %[output:348604a6] %[output:0e735326] %[output:32655eef] %[output:38a5c276] %[output:23eb441a] %[output:8aff8459] %[output:12af92de] %[output:5ab04c2d] %[output:6469fd3c] %[output:26559a5c] %[output:6156936f] %[output:706c5125] %[output:341926ac] %[output:31ab031b] %[output:6f34d664] %[output:545b86d3] %[output:5b6a8933] %[output:79f9df68] %[output:85b2da6e] %[output:89a00686] %[output:06e5ffcc] %[output:3019662b] %[output:6265b461] %[output:5e6d15eb] %[output:711ed926] %[output:176a1a49] %[output:6b4064a3] %[output:12abdea9] %[output:6f4d3815] %[output:8d49acd7] %[output:670a661f] %[output:5506aaa0] %[output:4f0de932] %[output:52c31584] %[output:72a8caef] %[output:1dfaddbd] %[output:3015739a] %[output:19230c7f] %[output:4ebac246] %[output:0782c66e] %[output:72dd4bab] %[output:50e4f630] %[output:9559d325] %[output:93d01694] %[output:77164893] %[output:27ec8cd7] %[output:4fd5fb83] %[output:2b10ba8d] %[output:4b2c9e37] %[output:5ee33073] %[output:078c14be] %[output:4d3d7188] %[output:32874822] %[output:420d3b5d] %[output:067e16a2] %[output:29e37ae4] %[output:6be3e630] %[output:649dac2b] %[output:28c1c108] %[output:5caab718] %[output:12d08640] %[output:98c1c0cc] %[output:6de7f332] %[output:568a297a] %[output:052586da] %[output:3179190b] %[output:0498429a] %[output:2eb48a50] %[output:0db968ab] %[output:96bdb53e] %[output:183c372d] %[output:6fd57591] %[output:01cacfb3] %[output:36e9588b] %[output:78d880ab] %[output:3d0cb05a] %[output:85e5f7b2] %[output:5e838c69] %[output:40cda491] %[output:0a1d04ff] %[output:31de50d8] %[output:348d9cba] %[output:67944c16] %[output:4c88cd3b] %[output:44a36b77] %[output:5d756a61] %[output:4bc6d90e] %[output:48461f64] %[output:8707e8dc] %[output:7dec6592] %[output:809c06c9] %[output:0d69a77e] %[output:62e51676] %[output:1e79beea] %[output:68dcd6a2] %[output:0e139b41] %[output:972cb268] %[output:7dee650a] %[output:53265cb6] %[output:2f836a3c] %[output:1279f925] %[output:1e79abd1] %[output:0871c531] %[output:711de2f7] %[output:849fd9f2] %[output:25187322] %[output:7a2db19e] %[output:1f449b17] %[output:470f987c] %[output:2dbf6996] %[output:8f5924d2] %[output:5fe95c90] %[output:61bc45ec] %[output:0200f8c6] %[output:4b9efe65] %[output:6fe38be9] %[output:671a4fe6] %[output:8490d017] %[output:453716cd] %[output:01377b30] %[output:057efccc] %[output:9d7b33a5] %[output:0a7e184a] %[output:5fe0a26d] %[output:2288493f] %[output:84cd1fcb] %[output:107cc1a2] %[output:24c5a872] %[output:86abaebb] %[output:88c38530] %[output:7eb7849b] %[output:335c9cd2] %[output:954b8e95] %[output:1559c752] %[output:2c754931] %[output:0089a340] %[output:8cf98f9f] %[output:9d88e280] %[output:310fb030] %[output:84c9d427] %[output:2c7ccea4] %[output:4ee8d51e] %[output:6701619c] %[output:19389fa5] %[output:7b020fc1] %[output:5cf875aa] %[output:2c0071ac] %[output:0e117794] %[output:1a108d4a] %[output:4941d0f0] %[output:27332270] %[output:83a38437] %[output:87ba3e0e] %[output:9354a170] %[output:8935b9a8] %[output:14827de9] %[output:4522cd14] %[output:5e47c4fd] %[output:1841baca] %[output:56f0f189] %[output:7d41a89c] %[output:4f11a50f] %[output:6118c065] %[output:1edf1d9d] %[output:92b757ba] %[output:007152c3] %[output:5130378c] %[output:4460b27f] %[output:6e71e5a6] %[output:9a8512ca] %[output:72b86ac9] %[output:5f2df76d] %[output:3cb32f7d] %[output:5619ac77] %[output:76051958] %[output:6852bfda] %[output:15ae0e44] %[output:20382e81] %[output:40e51824] %[output:225bfc75] %[output:87eda5ff] %[output:785cacc8] %[output:9c6c4a37] %[output:0f2f217a] %[output:94c8c355] %[output:789602f7] %[output:53c66c76] %[output:4590b788] %[output:07398029] %[output:0dab80da] %[output:9c748129] %[output:4bbc1d6a] %[output:260e5fa3] %[output:15ebee78] %[output:1367f547] %[output:89101d1c] %[output:47629766] %[output:435bb8c3] %[output:59c5305d] %[output:9fdb9ee0] %[output:1d77369a] %[output:817ef642] %[output:8c3bcdb0] %[output:9d91c16d] %[output:75b224ef] %[output:95fe3e0a] %[output:229a197c] %[output:094e1dfb] %[output:93c897b1] %[output:50454003] %[output:85664dac] %[output:9f317986] %[output:4c717dce] %[output:71c0fefa] %[output:596ff493] %[output:35dcb14d] %[output:9f07c8e4] %[output:563c42b7] %[output:34ecf5de] %[output:8cfa271f] %[output:6e1f2b22] %[output:77177098] %[output:1cc2ec99] %[output:09ca66c4] %[output:9e47b197] %[output:5d1edf3e] %[output:0c2a209f] %[output:84093a6f] %[output:71baf1eb] %[output:7cc41676] %[output:2e190aa3] %[output:31b8fb47] %[output:56946194] %[output:66704650] %[output:6fe207b8] %[output:3412f8d0] %[output:10d3d003] %[output:176ce60f] %[output:100c4b9d] %[output:69f0d83d] %[output:02f02f99] %[output:65a904d7] %[output:595487bb] %[output:33d77f95] %[output:0f233cb7] %[output:72b81796] %[output:1efed5b3] %[output:4afba3ed] %[output:3e1ea0ef] %[output:3f81ca1d] %[output:0de7c66e] %[output:94ea07b3] %[output:1263c7ee] %[output:3fb0e022] %[output:07548101] %[output:2bd6660b] %[output:9631a3b6] %[output:0fe32050] %[output:7e6dc31d] %[output:57fdea43] %[output:57eb0124] %[output:5c1c7eed] %[output:56da561a] %[output:49f6aaff] %[output:2d94ab1b] %[output:16ed10c7] %[output:647f8331] %[output:31f4acc8] %[output:71467855] %[output:78af054d] %[output:2e39546d] %[output:15c6d799] %[output:988c108e] %[output:80ad93c2] %[output:6541fe03] %[output:5800b3bd] %[output:3d35514e] %[output:630c9768] %[output:6696a670] %[output:2a89fa5d] %[output:2d7bae09] %[output:4c0c01ba] %[output:8b7a9749] %[output:82425e09] %[output:8c1ebc43] %[output:6e92b916] %[output:458c0924] %[output:35873879] %[output:8ecab688] %[output:1c815e97] %[output:2e48fc2b] %[output:0f546a4a] %[output:238a16a4] %[output:74055424] %[output:269d765f] %[output:66476c8e] %[output:8c09256b] %[output:4965f7d7] %[output:7e4b556b] %[output:1993e5ce] %[output:0b9215e5] %[output:4851e2bf] %[output:333af62d] %[output:74785caa] %[output:6ae14656] %[output:03afbbac] %[output:1f1522b0] %[output:8c6171b5] %[output:9a426651] %[output:191f5dfb] %[output:37066068] %[output:3a5b5539] %[output:4dd83f06] %[output:24eecd03] %[output:676b2c9e] %[output:3d659633] %[output:6b656174] %[output:6b6ae280] %[output:5a3ba033] %[output:2fe86fda] %[output:90c77023] %[output:2662aed1] %[output:0b1b1a6e] %[output:292235d2] %[output:1b060e02] %[output:4b08f57d] %[output:026d3b12] %[output:9b375755] %[output:604167cb] %[output:43cb3863] %[output:6c799ba7] %[output:31c088b3]

        labelMap = containers.Map('KeyType','double','ValueType','any');
        labelMap(rec.m_in) = 'm_in';
        labelMap(rec.m_Df) = 'm_Df';
        labelMap(rec.m_d)  = 'm_d';
        labelMap(rec.m_Ds) = 'm_Ds';
        labelMap(rec.m_c)  = 'm_c';
        labelMap(rec.m_Cs) = 'm_Cs';
        if hasMerge && ~isnan(rec.m_Nc)
            labelMap(rec.m_Nc) = 'm_Nc';
        end

        printRxnEquationTagged(model, rIdx, labelMap); %[output:1e2a13fd] %[output:5fc9c272] %[output:8e909393] %[output:6ce76db4] %[output:24f0d408] %[output:9e3a21a8] %[output:78aa188f] %[output:941ae897] %[output:7cd4819f] %[output:5d38cde1] %[output:96177af1] %[output:21daa189] %[output:2df32718] %[output:27985249] %[output:0238d820] %[output:7e54506d] %[output:46f47246] %[output:72943350] %[output:4fb5dec0] %[output:8383e761] %[output:916b74be] %[output:670b7101] %[output:90263e34] %[output:4d48e02b] %[output:27e181e8] %[output:1bb8047b] %[output:531d9a8f] %[output:387c81e3] %[output:48a815cc] %[output:260f8ed6] %[output:42e8b30c] %[output:5d80f93a] %[output:24691816] %[output:9512e15a] %[output:11488d6f] %[output:828cceab] %[output:1458876e] %[output:670129ba] %[output:03b3519e] %[output:458d9f55] %[output:0d417563] %[output:1313aace] %[output:89a67045] %[output:4a0736d7] %[output:7771c822] %[output:88b89ccb] %[output:85f37d9e] %[output:88126fb7] %[output:19daaab7] %[output:4f55ef65] %[output:9089a3df] %[output:130d0db9] %[output:616634f8] %[output:5e16a1a6] %[output:1f4d44d5] %[output:195e79e9] %[output:8d471e8a] %[output:59ea689e] %[output:9b4df873] %[output:2020556b] %[output:71edfd8a] %[output:0b61021a] %[output:6a64be39] %[output:60476223] %[output:46c680c0] %[output:0a3482bc] %[output:4b509dc2] %[output:827cef35] %[output:69f0ad63] %[output:00821ae3] %[output:951e3170] %[output:20a17029] %[output:17bb76d4] %[output:3f41f5fc] %[output:3437bdcc] %[output:84249311] %[output:137001c3] %[output:0607afd0] %[output:167ba8ec] %[output:0524fd24] %[output:8cfe5244] %[output:652982c3] %[output:975edb19] %[output:9a36303a] %[output:7d899b6d] %[output:70cfb4b8] %[output:8b59263f] %[output:9c2ef04e] %[output:5b4fedb3] %[output:26c34093] %[output:78f35f1a] %[output:0741cc74] %[output:2abf5672] %[output:94cd7456] %[output:6db8f295] %[output:0175c56b] %[output:472a64f7] %[output:5eaf5562] %[output:271e9aa0] %[output:49b57f06] %[output:4004b77c] %[output:8a3b1f01] %[output:3a46bd0a] %[output:0043c767] %[output:5381724c] %[output:34f29fbc] %[output:5fecbd07] %[output:23048000] %[output:28583c2c] %[output:24dfe9ca] %[output:30085f2c] %[output:060e2f3b] %[output:613cb2ac] %[output:18d4a4f6] %[output:8c5a63e6] %[output:4acad47d] %[output:4d224d46] %[output:6c3f6104] %[output:72c81ce2] %[output:18158001] %[output:2ac7d9db] %[output:40150084] %[output:18e69d18] %[output:89eaac4b] %[output:66a1b145] %[output:3baedd67] %[output:6e78db42] %[output:832cd459] %[output:7c73bf20] %[output:26bb3397] %[output:9583328a] %[output:919998b1] %[output:703c4507] %[output:28b8e5df] %[output:84c73fa4] %[output:15e108a3] %[output:0871038a] %[output:07198641] %[output:80227d2a] %[output:78747869] %[output:34bcf3b9] %[output:8bbb5994] %[output:7a1b5f28] %[output:4a408b36] %[output:11aed6ee] %[output:0470a6f3] %[output:5afbe187] %[output:721e3de0] %[output:83380371] %[output:982f265d] %[output:3ea0802b] %[output:929e1ca3] %[output:43b88841] %[output:1ffb7800] %[output:60e04d4b] %[output:44324e49] %[output:31e83f45] %[output:477bdec0] %[output:9b981943] %[output:78e0d8cb] %[output:4ca37162] %[output:71ab5a36] %[output:5879e098] %[output:4d930ac9] %[output:45ef46e5] %[output:35269413] %[output:05d3d0cc] %[output:2a2f66e1] %[output:4870c334] %[output:0a613f40] %[output:8d4725f1] %[output:5622632e] %[output:71a83d36] %[output:5954ce2b] %[output:257d2b8a] %[output:9ca9b3b1] %[output:017b4ce0] %[output:3c18bc7a] %[output:0a731811] %[output:2b86b152] %[output:9c723b70] %[output:8144d6c2] %[output:13f86f3e] %[output:9629fa64] %[output:7d25b0e0] %[output:5b45c4af] %[output:77245c68] %[output:874b7acd] %[output:6a607bda] %[output:8b87cd21] %[output:8cb335a2] %[output:449791b6] %[output:3d541ac5] %[output:693746f5] %[output:74f69307] %[output:3bdadd49] %[output:615f906d] %[output:853afc98] %[output:1b4a6c9e] %[output:77f4d03a] %[output:90999d2f] %[output:29e21551] %[output:1fa89868] %[output:9a3a39e5] %[output:35939cee] %[output:595cade5] %[output:1ffefa54] %[output:01bcadbe] %[output:8c150137] %[output:4d74db61] %[output:1b67ca91] %[output:553f6c3a] %[output:3403f4bb] %[output:6c65914e] %[output:0ccdf30e] %[output:979e90c9] %[output:4a0039cd] %[output:34492b08] %[output:6b95b09f] %[output:9f596d47] %[output:100ab65d] %[output:494a66fe] %[output:631c4891] %[output:3931a7fe] %[output:45c6b1cb] %[output:3bcc992f] %[output:084ae19c] %[output:0e7bf32f] %[output:77391a3d] %[output:5d80bbf9] %[output:718dfdf0] %[output:3107230c] %[output:1c7d0743] %[output:756b2e28] %[output:066cb8ca] %[output:639977bf] %[output:9c09908d] %[output:198ba20e] %[output:2983c80a] %[output:264c9950] %[output:94d31078] %[output:8ac45805] %[output:22459553] %[output:98247548] %[output:52ee7670] %[output:33300347] %[output:4835e076] %[output:3a52db53] %[output:88d49548] %[output:1eb13aa2] %[output:55512952] %[output:55e57b72] %[output:6c65d7d7] %[output:6285a20c] %[output:45103c3f] %[output:363dc789] %[output:880eba82] %[output:5b385c24] %[output:5a0b3f59] %[output:3eb6822a] %[output:849c2831] %[output:939fd7d4] %[output:06a56b46] %[output:2a1a733b] %[output:547c88d6] %[output:166daa85] %[output:29d6354f] %[output:915d5299] %[output:20160692] %[output:582ea795] %[output:4ad23462] %[output:42735cf0] %[output:7bedcd25] %[output:894d3deb] %[output:10286563] %[output:7542eeac] %[output:3a22683e] %[output:0723ccaa] %[output:8dd7d69d] %[output:3dce8aa4] %[output:00d4f58c] %[output:4988c375] %[output:81540859] %[output:3722b70b] %[output:7337673a] %[output:9982a720] %[output:93014233] %[output:904c8642] %[output:0d758d2d] %[output:978cbff1] %[output:25ddbc36] %[output:15caa544] %[output:44e4dd66] %[output:9519a2c5] %[output:8c986a45] %[output:417226df] %[output:8acacfdf] %[output:726d575b] %[output:724a0098] %[output:087d6db6] %[output:98aed323] %[output:27209e67] %[output:9912d86a] %[output:2be5abb0] %[output:9b49b12a] %[output:7e59b6c1] %[output:122cc42c] %[output:6abad913] %[output:9cb00fba] %[output:2780cd19] %[output:9afbdb24] %[output:9328dd14] %[output:202c22cf] %[output:10ca7a85] %[output:6f6bd5ec] %[output:9d260a07] %[output:7f972642] %[output:5fbfd73e] %[output:2a223d6d] %[output:4f97d10b] %[output:6b733f1f] %[output:6c2758f6] %[output:13093f51] %[output:87df82db] %[output:15df4db2] %[output:3d10e408] %[output:3fce5783] %[output:83b85de5] %[output:3f1778a2] %[output:3df9d41e] %[output:28242d5a] %[output:2eced40b] %[output:998cfb83] %[output:2df77f7d] %[output:621cb511] %[output:7eec5c98] %[output:5a754beb] %[output:6b15f1bf] %[output:188bb068] %[output:4697bd6a] %[output:3a27bc2d] %[output:32cba3d9] %[output:92a08e0d] %[output:0dead9ec] %[output:00827daf] %[output:1681755e] %[output:1d5a3c80] %[output:92753874] %[output:5f7cf004] %[output:2c68e2fe] %[output:63d3eca3] %[output:0724f1fd] %[output:66cadb2b] %[output:5bb5f5dd] %[output:8d4004b5] %[output:1b27e458] %[output:37a19b55] %[output:51e41887] %[output:2bb54fae] %[output:1979521d] %[output:55ab21fc] %[output:8ee465a8] %[output:4db6c03f] %[output:9d5ac708] %[output:31ff35db] %[output:34fd466c] %[output:325a5cea] %[output:4442ac63] %[output:2398d181] %[output:3b25e4ed] %[output:353cefec] %[output:1cacf52c] %[output:7f70873e] %[output:8365cfde] %[output:4d97eb59] %[output:204c64f6] %[output:7fff8512] %[output:5663a411] %[output:2cb0aad9] %[output:5879c3d1] %[output:9cb315a1] %[output:4affa2b7] %[output:40639224] %[output:4fd9c140] %[output:106c3f5f] %[output:5fb140b2] %[output:5758c6dc] %[output:4f501cef] %[output:7cf86e65] %[output:0a9bd84c] %[output:523828f5] %[output:98df5192] %[output:39b2c8bc] %[output:17a677a2] %[output:5b36837b] %[output:43695f90] %[output:59eac599] %[output:8f612036] %[output:0634f861] %[output:1c596e5a] %[output:16615950] %[output:5e5b0577] %[output:45a23646] %[output:576ced0c] %[output:6ea69738] %[output:4b4f8b16] %[output:92c70104] %[output:103f13a7] %[output:7a288143] %[output:92caddf3] %[output:31cf8627] %[output:92895d9e] %[output:802f09e7] %[output:7cdf7b9c] %[output:371a5c69] %[output:0741df89] %[output:861d6497] %[output:2c69ed35] %[output:6bb095b0] %[output:68af04fa] %[output:1bd58581] %[output:23a9c9b5] %[output:88a3d7fd] %[output:1d084f55] %[output:07b88109] %[output:28427410] %[output:665a34aa] %[output:106a5974] %[output:1d02c128] %[output:89b1df8d] %[output:46dc8241] %[output:8a908e03] %[output:50d063fb] %[output:146eabd3] %[output:992cfda3] %[output:718d5e32] %[output:31e5a39c] %[output:573c252c] %[output:168f2244] %[output:72871d37] %[output:10e7ee7d] %[output:60f3dea8] %[output:5db57258] %[output:9ea0c718] %[output:066cc06c] %[output:1755fa69] %[output:24505d46] %[output:2455ce65] %[output:68c6a51a] %[output:37a6723b] %[output:0571f99e] %[output:2edab454] %[output:986582c7] %[output:2f6836d5] %[output:50012c9a] %[output:55d7c6ed] %[output:2031a5ef] %[output:9590f45d] %[output:4a599c28] %[output:37951ece] %[output:64350e43] %[output:8652767b] %[output:4dfc6a23] %[output:7117daa6] %[output:13fe025f] %[output:545fbb0d] %[output:76b6bbc7] %[output:8d72c36c] %[output:8e6eb953] %[output:56cdc211] %[output:427690b2] %[output:615610f2] %[output:49a2e1d1] %[output:0c4b77ac] %[output:2cebc2b6] %[output:4d444e2e] %[output:36f25790] %[output:62aaa7f7] %[output:712c7708] %[output:223239dd]

    end
end %[output:group:30832ef2]

if nF == 0
    fprintf('\nNo joint diff+counter+NN motifs available to print ');
    fprintf('(check that nnStep1 / nnStep3 are non-empty).\n');
end

if nF > 0 %[output:group:58d88595]
    FullMotifTable = struct2table(fullJoinRows);
    save('nn_full_motifs.mat', 'FullMotifTable', 'fullJoinRows', 'hasMerge');
    fprintf('\nSaved nn_full_motifs.mat (%d full motifs, ranked by gene complexity)\n', nF); %[output:6288ed4c]
end %[output:group:58d88595]



%%\


%%
%[text] ## ---------------- Helper functions ----------------
function ec = safeECCode(model, rIdx)
% Returns the EC code(s) for reaction rIdx, or "(none annotated)" if empty.
ec = "(none annotated)";
if isfield(model,'eccodes') && numel(model.eccodes) >= rIdx && strlength(strtrim(string(model.eccodes{rIdx}))) > 0
    ec = string(model.eccodes{rIdx});
end
end

function printRxnEquationTagged(model, rIdx, labelMap)
% Prints "coeff Formula(tag) + coeff Formula(tag) -> ..." for reaction
% rIdx, where (tag) is the motif notation (m_in, m_d, m_c, ...) for any
% metabolite present in labelMap. Metabolites not in labelMap print with
% no tag at all (just "coeff Formula"). No color, no external dependency.
S = model.S;
col = full(S(:,rIdx));
subIdx = find(col < 0);
prodIdx = find(col > 0);

fprintf('      Eqn: ');
printSideTagged(model, subIdx, col, true, labelMap);
fprintf('  ->  ');
printSideTagged(model, prodIdx, col, false, labelMap);
fprintf('\n');
end

function printSideTagged(model, idxList, col, isSubstrate, labelMap)
if isempty(idxList)
    fprintf('(none)');
    return;
end
for i = 1:numel(idxList)
    m = idxList(i);
    if isSubstrate, coeff = -col(m); else, coeff = col(m); end
    formula = safeMetFormula(model, m);

    if abs(coeff - 1) < 1e-6
        coeffStr = '';
    else
        coeffStr = sprintf('%.3g ', coeff);
    end

    if isKey(labelMap, m)
        token = sprintf('%s%s(%s)', coeffStr, formula, labelMap(m));
    else
        token = sprintf('%s%s', coeffStr, formula);
    end
    fprintf('%s', token);
    if i < numel(idxList)
        fprintf(' + ');
    end
end
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

function f = safeMetFormula(model, mIdx)
% Returns the metabolite formula if the field exists, else empty string.
f = "";
if isfield(model,'metFormulas') && numel(model.metFormulas) >= mIdx && ~isempty(model.metFormulas{mIdx})
    f = string(model.metFormulas{mIdx});
end
end
 
function g = safeGPR(model, rIdx)
% Returns the GPR rule string, or "(none annotated)" if empty.
g = "(none annotated)";
if isfield(model,'grRules') && numel(model.grRules) >= rIdx && strlength(strtrim(string(model.grRules{rIdx}))) > 0
    g = string(model.grRules{rIdx});
end
end
 
function eqStr = safeRxnEquation(model, rIdx)
% Builds a simple "substrates -> products" string using metabolite IDs
% and stoichiometric coefficients, reading directly from model.S. This
% avoids depending on any external helper function (e.g.
% reactionToString_chm) that may or may not be on the path.
S = model.S;
col = full(S(:,rIdx));
subIdx = find(col < 0);
prodIdx = find(col > 0);
 
subStr = strjoin(arrayfun(@(m) sprintf('%.3g %s', -col(m), safeMetFormula(model, m)), ...
                           subIdx, 'UniformOutput', false), ' + ');
prodStr = strjoin(arrayfun(@(m) sprintf('%.3g %s', col(m), safeMetFormula(model, m)), ...
                            prodIdx, 'UniformOutput', false), ' + ');
if strlength(subStr) == 0, subStr = "(none)"; end
if strlength(prodStr) == 0, prodStr = "(none)"; end
eqStr = sprintf('%s  ->  %s', subStr, prodStr);
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":52}
%---
%[output:986f0ab3]
%   data: {"dataType":"text","outputData":{"text":"Adding matlab path to: \/home\/ccs\/torres-gomez\/Downloads\/bioelectrodynamics\/code\n","truncated":false}}
%---
%[output:6d612b06]
%   data: {"dataType":"text","outputData":{"text":"\n\n      _____   _____   _____   _____     _____     |\n     \/  ___| \/  _  \\ |  _  \\ |  _  \\   \/ ___ \\    |   COnstraint-Based Reconstruction and Analysis\n     | |     | | | | | |_| | | |_| |  | |___| |   |   The COBRA Toolbox - 2026\n     | |     | | | | |  _  { |  _  \/  |  ___  |   |\n     | |___  | |_| | | |_| | | | \\ \\  | |   | |   |   Documentation:\n     \\_____| \\_____\/ |_____\/ |_|  \\_\\ |_|   |_|   |   <a href=\"http:\/\/opencobra.github.io\/cobratoolbox\">http:\/\/opencobra.github.io\/cobratoolbox<\/a>\n                                                  | \n\n > Checking if git is installed ...  Done (version: 2.30.2).\n > Checking if the repository is tracked using git ...  Done.\n > Checking if curl is installed ...  Done.\n > Checking if remote can be reached ...  Done.\n > Initializing and updating submodules (this may take a while)... Done.\n > Adding all the files of The COBRA Toolbox ...  Done.\n > Define CB map output... set to svg.\n > TranslateSBML is installed and working properly.\n > Configuring solver environment variables ...\n   - [*---] ILOG_CPLEX_PATH: --> set this path manually after installing the solver ( see <a href=\"https:\/\/opencobra.github.io\/cobratoolbox\/docs\/solvers.html\">instructions<\/a> )\n   - [*---] GUROBI_PATH: --> set this path manually after installing the solver ( see <a href=\"https:\/\/opencobra.github.io\/cobratoolbox\/docs\/solvers.html\">instructions<\/a> )\n   - [*---] TOMLAB_PATH: --> set this path manually after installing the solver ( see <a href=\"https:\/\/opencobra.github.io\/cobratoolbox\/docs\/solvers.html\">instructions<\/a> )\n   - [*---] MOSEK_PATH: --> set this path manually after installing the solver ( see <a href=\"https:\/\/opencobra.github.io\/cobratoolbox\/docs\/solvers.html\">instructions<\/a> )\n   Done.\n > Checking available solvers and solver interfaces ...Could not find installation of dqqMinos, so it cannot be tested\nGurobi installed at this location? \nLicence file current? \nCould not find installation of mosek, so it cannot be tested\nCould not find installation of quadMinos, so it cannot be tested\nCould not find installation of tomlab_snopt, so it cannot be tested\n Done.\n > Setting default solvers ...Could not find installation of mosek, so it cannot be tested\nCould not find installation of mosek, so it cannot be tested\n Done.\n > Saving the MATLAB path ... Done.\n   - The MATLAB path was saved as ~\/pathdef.m.\n\n > Summary of available solvers and solver interfaces\n\n\t\t\tSupport \t   LP \t MILP \t   QP \t MIQP \t  NLP \t   EP \t  CLP\n\t------------------------------------------------------------------------------\n\tdqqMinos     \tactive        \t    0 \t    - \t    0 \t    - \t    - \t    - \t    -\n\tglpk         \tactive        \t    1 \t    1 \t    - \t    - \t    - \t    - \t    -\n\tgurobi       \tactive        \t    0 \t    0 \t    0 \t    0 \t    - \t    - \t    -\n\tlp_solve     \tlegacy        \t    1 \t    - \t    - \t    - \t    - \t    - \t    -\n\tmatlab       \tactive        \t    1 \t    - \t    - \t    - \t    1 \t    - \t    -\n\tmosek        \tactive        \t    0 \t    - \t    0 \t    - \t    - \t    0 \t    0\n\tpdco         \tactive        \t    1 \t    - \t    1 \t    - \t    - \t    1 \t    -\n\tqpng         \tpassive       \t    - \t    - \t    1 \t    - \t    - \t    - \t    -\n\tquadMinos    \tactive        \t    0 \t    - \t    - \t    - \t    - \t    - \t    -\n\ttomlab_snopt \tpassive       \t    - \t    - \t    - \t    - \t    0 \t    - \t    -\n\t------------------------------------------------------------------------------\n\tTotal        \t-             \t    4 \t    1 \t    2 \t    0 \t    1 \t    1 \t    0\n\n + Legend: - = not applicable, 0 = solver not compatible or not installed, 1 = solver installed.\n\n\n > You can solve LP problems using: 'glpk' - 'pdco' \n > You can solve MILP problems using: 'glpk' \n > You can solve QP problems using: 'pdco' \n > You can solve MIQP problems using: \n > You can solve NLP problems using: \n > You can solve EP problems using: 'pdco' \n > You can solve CLP problems using: \n\nGurobi installed at this location? \nLicence file current? \nGurobi installed at this location? \nLicence file current? \nGurobi installed at this location? \nLicence file current? \nGurobi installed at this location? \nLicence file current? \n> Checking for available updates ... skipped\nremoving: \/home\/ccs\/torres-gomez\/Downloads\/bioelectrodynamics\/code\/cobratoolbox~\/src\/analysis\/thermo\/componentContribution\/new\nremoving: \/home\/ccs\/torres-gomez\/Downloads\/bioelectrodynamics\/code\/cobratoolbox~\/src\/analysis\/thermo\/groupContribution\/new\nremoving: \/home\/ccs\/torres-gomez\/Downloads\/bioelectrodynamics\/code\/cobratoolbox~\/src\/analysis\/thermo\/inchi\/new\nremoving: \/home\/ccs\/torres-gomez\/Downloads\/bioelectrodynamics\/code\/cobratoolbox~\/src\/analysis\/thermo\/molFiles\/new\nremoving: \/home\/ccs\/torres-gomez\/Downloads\/bioelectrodynamics\/code\/cobratoolbox~\/src\/analysis\/thermo\/protons\/new\nremoving: \/home\/ccs\/torres-gomez\/Downloads\/bioelectrodynamics\/code\/cobratoolbox~\/src\/analysis\/thermo\/trainingModel\/new\n","truncated":false}}
%---
%[output:0bb3a482]
%   data: {"dataType":"text","outputData":{"text":"\/home\/ccs\/torres-gomez\/Downloads\/bioelectrodynamics\/code\/cobratoolbox~\/initCobraToolbox.m\n","truncated":false}}
%---
%[output:16e16ec5]
%   data: {"dataType":"text","outputData":{"text":"Each model.subSystems{x} has been changed to a character array.\n","truncated":false}}
%---
%[output:6f744176]
%   data: {"dataType":"text","outputData":{"text":"Counter search on iMM904: 2806 mets, 4131 rxns.\n","truncated":false}}
%---
%[output:2b5d9fae]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"7×1 string array","name":"ans","rows":7,"type":"string","value":[["(2-amino-4-hydroxy-7,8-dihydropteridin-6-yl)methyl trihydrogen diphosphate"],["3-(imidazol-4-yl)-2-oxopropyl dihydrogen phosphate"],["H+"],["hydrogen cyanide"],["hydrogen peroxide"],["hydrogen sulfide"],["orotidine 5'-(dihydrogen phosphate)"]]}}
%---
%[output:937c7c33]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"3×1 string array","name":"ans","rows":3,"type":"string","value":[["5-phosphoribosyl-ATP"],["ATP"],["dATP"]]}}
%---
%[output:8934d6ad]
%   data: {"dataType":"text","outputData":{"text":"Enzyme-mediated rxns: 2709 out of 4131\n","truncated":false}}
%---
%[output:391715b2]
%   data: {"dataType":"text","outputData":{"text":"Total candidates (4-reaction topology): 0\n","truncated":false}}
%---
%[output:82491b5b]
%   data: {"dataType":"text","outputData":{"text":"\n\n############################################################\n","truncated":false}}
%---
%[output:193099ec]
%   data: {"dataType":"text","outputData":{"text":"FULL DIFFERENTIATOR + COUNTER + NN MOTIFS\n","truncated":false}}
%---
%[output:3a8cd22e]
%   data: {"dataType":"text","outputData":{"text":"Showing 5 simplest of 48 scored (ranked by total distinct genes, ascending)\n","truncated":false}}
%---
%[output:2be4d4c6]
%   data: {"dataType":"text","outputData":{"text":"############################################################\n","truncated":false}}
%---
%[output:90480593]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:2f5ab102]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 1  (diff_idx=245, cnt_idx=4, total genes=11)\n","truncated":false}}
%---
%[output:31d9fd73]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:2a57d334]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:6fd3818a]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:99e6fccd]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:964fb803]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:579df49a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:04de2ea6]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:1e2a13fd]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:87d38746]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:72e5a225]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:92c5e06b]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:5fc9c272]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:962adbdf]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:6e54ee94]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:2fa6c30a]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:8e909393]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:1e0be3a5]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:65a81a21]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:8b427032]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:6ce76db4]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:59b39ee9]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:794d8672]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:4ba27199]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:24f0d408]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:30f0434f]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:5c949d59]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:925d121e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:9e3a21a8]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:2ad6032f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:6b91d1ff]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:981320ab]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:78aa188f]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:73591e18]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:228a90b3]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:6403bf32]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:941ae897]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:441b3eb8]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:762ff6e3]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:5d753d9b]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:7cd4819f]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:34d2911b]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:1ae71e83]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:8bb958de]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:5d38cde1]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C11H15N5O3S + HO4P(m_Nc)  ->  C5H5N5 + C6H11O7PS\n","truncated":false}}
%---
%[output:6700baa7]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:3aefed9b]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:20d37575]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:5a3d340d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:30d2eb08]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:9b29afdd]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:9014db65]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:50c8e3c7]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:2029ab25]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.28\n","truncated":false}}
%---
%[output:96177af1]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:4dc2e804]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:3802177c]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:3e2821ca]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:21daa189]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:93d5880c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:66364b24]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:69ce1586]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:2df32718]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:12c73b99]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:43bdd8ca]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:6b0a14bd]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR017W\n","truncated":false}}
%---
%[output:27985249]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:2f7409df]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:2f733ce7]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:7d4c6c8e]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:0238d820]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:90981ff5]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:5408d0af]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:9a784dc8]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:7e54506d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:86460037]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:148ad159]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:99b2f8a2]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:46f47246]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:5fb059d9]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:6be7e79a]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:965bbc39]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0075         | 5'-methylthioadenosine phosphorylase\n","truncated":false}}
%---
%[output:72943350]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:2645637c]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:2d9146c8]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:44807ecc]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:4fb5dec0]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:0c743432]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:46002d9c]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:7b619fd0]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:8383e761]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H6O27P7 + 2 H + HO4P(m_Nc)  ->  C6H7O30P8 + H2O\n","truncated":false}}
%---
%[output:808bd563]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:74fa3691]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:9209d04c]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:4c628799]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:22a8bfc7]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:22a201bf]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:6d83ae88]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:3dac186d]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:8f6b2692]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:916b74be]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:7fd78e3f]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:93ab954d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:849175bc]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:670b7101]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:657192d9]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:578ffdb6]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:5007c0d9]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:90263e34]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:3c13c830]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:47f98be8]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:722d9ba6]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:4d48e02b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:4be7c015]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.7.4.21;2.7.4.24\n","truncated":false}}
%---
%[output:9c1576bd]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:99f6fe94]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:27e181e8]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:9fa608e5]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:1c3157c9]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:94e2264f]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:1bb8047b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:62a42c84]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:70bbf2e1]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:77b7c0f3]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:531d9a8f]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:0e1e59fe]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:0de14ea5]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:3b45ba51]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:387c81e3]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:2af384bc]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:79f2e308]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:51f71ce4]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:48a815cc]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:03c1557f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:619957f0]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:92b262fc]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:260f8ed6]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H6O27P7 + 2 H + HO4P(m_Nc)  ->  C6H7O30P8 + H2O\n","truncated":false}}
%---
%[output:8f1ada4a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:8fcbaeb6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:7f48dbcb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:9733901e]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:3749cea5]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:687bcc2b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR410W\n","truncated":false}}
%---
%[output:8cc72d13]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 2  (diff_idx=245, cnt_idx=4, total genes=11)\n","truncated":false}}
%---
%[output:1acc2978]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:24e4049a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:42e8b30c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:2bbb1da1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:5512b39e]
%   data: {"dataType":"text","outputData":{"text":"\n============================================================\n","truncated":false}}
%---
%[output:56951bb7]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 3  (diff_idx=245, cnt_idx=4, total genes=11)\n","truncated":false}}
%---
%[output:5d80f93a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:0ee78ee8]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:5c73fc0a]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 4  (diff_idx=245, cnt_idx=4, total genes=11)\n","truncated":false}}
%---
%[output:91575452]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:24691816]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:1cabc945]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 5  (diff_idx=245, cnt_idx=4, total genes=11)\n","truncated":false}}
%---
%[output:2a881c51]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:7afdbf40]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 6  (diff_idx=245, cnt_idx=4, total genes=11)\n","truncated":false}}
%---
%[output:9512e15a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:9dbc1e29]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:96d3faf8]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 7  (diff_idx=245, cnt_idx=4, total genes=11)\n","truncated":false}}
%---
%[output:75f2adc8]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:11488d6f]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:2f9e96e5]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:9402ba8b]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 8  (diff_idx=245, cnt_idx=4, total genes=11)\n","truncated":false}}
%---
%[output:8b37fa4e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:828cceab]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:3fac79ef]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:49026482]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 9  (diff_idx=245, cnt_idx=4, total genes=11)\n","truncated":false}}
%---
%[output:6d41f771]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 10  (diff_idx=245, cnt_idx=4, total genes=11)\n","truncated":false}}
%---
%[output:1458876e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:619b6231]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:04458a67]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.7.4.21;2.7.4.24\n","truncated":false}}
%---
%[output:5a9b1eb4]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:670129ba]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:25b6f40b]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 11  (diff_idx=245, cnt_idx=4, total genes=11)\n","truncated":false}}
%---
%[output:1165b0b4]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 12  (diff_idx=245, cnt_idx=4, total genes=11)\n","truncated":false}}
%---
%[output:96147fd6]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 13  (diff_idx=245, cnt_idx=4, total genes=11)\n","truncated":false}}
%---
%[output:03b3519e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:852cd886]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:06319b5f]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0088         | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:14b2fe30]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 14  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:458d9f55]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + 2 H + C6H7O21P5 + 2 HO4P(m_Nc)  ->  C6H7O30P8 + C10H12N5O10P2 + 2 H2O\n","truncated":false}}
%---
%[output:6ec3eea5]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 15  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:7ca31de2]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:6b36469d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:6a7010e6]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 16  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:60a6c84c]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 17  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:4b1d3df3]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 18  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:08f05f0e]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:9e5ebd64]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:465b8e3a]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 19  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:0d417563]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:6c5b0824]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 20  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:2b19d544]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:9ef7df01]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 21  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:1313aace]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:972a5e47]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:50a650c0]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 22  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:08a1a377]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:89a67045]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:8b582948]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:63cafd0e]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 23  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:8db61185]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:4a0736d7]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:349b76d6]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:6ecd9066]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 24  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:5023a33e]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 25  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:7771c822]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:7681edb9]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:6469a447]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:43084eaf]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:88b89ccb]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:971bbb8d]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 26  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:3725d90c]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 27  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:681ba772]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:85f37d9e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:5d4722e5]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 28  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:4b002ffc]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR410W\n","truncated":false}}
%---
%[output:4c48355e]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:88126fb7]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:13590dff]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:819564cf]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 29  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:7dd8c16a]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 30  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:19daaab7]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:57561bf6]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 31  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:0dd8ed02]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:64c521c2]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:4f55ef65]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H10O5 + HO4P(m_Nc)  ->  C6H11O9P\n","truncated":false}}
%---
%[output:4227a215]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 32  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:773f39b2]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 33  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:9c84642a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.7.4.21\n","truncated":false}}
%---
%[output:4573bb35]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 34  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:5a005d6a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:3e6b709b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:732834f8]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 35  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:91c490cb]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 36  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:236518b4]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 37  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:9089a3df]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:1dad3065]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:1ce80a59]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:5d71f650]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 38  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:130d0db9]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:6450ea47]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:10685e5c]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 39  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:8df2beb6]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:616634f8]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:02e5eb98]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 40  (diff_idx=245, cnt_idx=4, total genes=12)\n","truncated":false}}
%---
%[output:388b3e5d]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:648b807c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:5e16a1a6]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:574176bd]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 41  (diff_idx=245, cnt_idx=4, total genes=13)\n","truncated":false}}
%---
%[output:9ae86943]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:7d35be0a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:1f4d44d5]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:50adb809]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:883afde2]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 42  (diff_idx=245, cnt_idx=4, total genes=13)\n","truncated":false}}
%---
%[output:95583f69]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:195e79e9]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:6c4422c6]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:97dbf21c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:1cf4ac40]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:8d471e8a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:1f4d246d]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 43  (diff_idx=245, cnt_idx=4, total genes=13)\n","truncated":false}}
%---
%[output:1ae29f36]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:872455c6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:59ea689e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:46b4e7e3]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:0af20fff]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:3756a6f6]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:9b4df873]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:995a343b]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 44  (diff_idx=245, cnt_idx=4, total genes=14)\n","truncated":false}}
%---
%[output:45ff7400]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:495b55aa]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:2020556b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C11H15N2O5 + HO4P(m_Nc)  ->  C5H9O8P + H + C6H6N2O\n","truncated":false}}
%---
%[output:0c7bb582]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:657a287f]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 45  (diff_idx=245, cnt_idx=4, total genes=14)\n","truncated":false}}
%---
%[output:51a8354f]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0089         | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:72f404a4]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 46  (diff_idx=228, cnt_idx=5, total genes=14)\n","truncated":false}}
%---
%[output:8bc95633]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:9786816e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:79224839]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 47  (diff_idx=228, cnt_idx=5, total genes=14)\n","truncated":false}}
%---
%[output:469aff16]
%   data: {"dataType":"text","outputData":{"text":" MOTIF 48  (diff_idx=677, cnt_idx=1, total genes=15)\n","truncated":false}}
%---
%[output:2f9043e6]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:71edfd8a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:96d34659]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:92aa1e38]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:14b2ac8e]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:0b61021a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:13ad370a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.1.1\n","truncated":false}}
%---
%[output:4bfebf99]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:51fba45e]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:6a64be39]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:03440e61]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YDR017C\n","truncated":false}}
%---
%[output:346dd419]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:057708ae]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:60476223]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:6cefafbd]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:0a8a56ae]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:23c8360f]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:46c680c0]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:6d6749b8]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:58b1fce7]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:654dbb6a]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:0a3482bc]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:3897453e]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:1e78f4da]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:6f362ef1]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:4b509dc2]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:64f56ee7]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:48868e7b]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:0ea0a0a7]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:827cef35]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:9dd8af12]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:8008f20d]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:4d8a3c9b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:69f0ad63]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:00f8fdba]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:011b3745]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:450e25c1]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:00821ae3]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H13N5O5 + HO4P(m_Nc)  ->  C5H9O8P + C5H5N5O\n","truncated":false}}
%---
%[output:24a88e5f]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:9972d2fa]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:0a0c449f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:5d5596df]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:33cdc236]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:1ab575f4]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:64f71a5f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:4125fb39]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:620da475]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:951e3170]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:48b34426]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:987b6de9]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:51c89d90]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:20a17029]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:87594fe4]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:82bc4190]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:33d0cf16]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:17bb76d4]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:6fa8478e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:9e286bcd]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:38c0842e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:3f41f5fc]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:4994efd2]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:124e295c]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:18e27761]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:3437bdcc]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:5075e836]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:5b84f492]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.1;2.4.2.28;3.2.2.3\n","truncated":false}}
%---
%[output:98f5e7f9]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:84249311]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:081a0053]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:86689416]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:63cb3ce2]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:137001c3]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:0961979e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:900ccf5f]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:2bfd9fc8]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:0607afd0]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:0682ab7b]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:806f2c93]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:88d7adb4]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:167ba8ec]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:61b974e3]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:818abd41]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:45b123e5]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:0524fd24]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N4O5 + HO4P(m_Nc)  ->  C5H9O8P + C5H4N4O\n","truncated":false}}
%---
%[output:047ba959]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:96779ef6]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR160W\n","truncated":false}}
%---
%[output:279671dc]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:20f8b587]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:24f7dd0f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:837c2821]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0358         | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:6175d1a8]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:17edb768]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:9bb19f9b]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:8cfe5244]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:2214a907]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:8356dff8]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:5af7bd9c]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:652982c3]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:81724c1a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:805f5f21]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:23071a4f]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:975edb19]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:461543db]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:6cdbdd1e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:018c4b56]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:9a36303a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:23060af2]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:00dac4ee]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:9c4e81d9]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:7d899b6d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:8dd47c55]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:8bc7e3e3]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:0a514dff]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:70cfb4b8]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:2ca73f8e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:19163115]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:830360bd]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:8b59263f]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:61c25090]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:37a33a25]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:2d764a1e]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:9c2ef04e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:892cba43]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:94afa365]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:5a4088d1]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.1\n","truncated":false}}
%---
%[output:5b4fedb3]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:73cabbd7]
%   data: {"dataType":"text","outputData":{"text":"============================================================\n","truncated":false}}
%---
%[output:109ff7bc]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:6bc6e900]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:26c34093]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O10PS + HO4P(m_Nc)  ->  C10H12N5O10P2 + H + O4S\n","truncated":false}}
%---
%[output:694e7e75]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:948d54b0]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:08bf6b4d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:3067d587]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:1155f6ef]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:6c0c4b6b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:80c159c0]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:4bcd9333]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:59987011]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:78f35f1a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:6cf04fc7]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:0d79c1ff]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:67cfd80a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:0741cc74]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:6cb8ddac]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:7c442532]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:7cbf001a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:2abf5672]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:767e01dc]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:640e8294]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:7a190c11]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:94cd7456]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:316a5dbe]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:0f32864a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:2c44f614]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:6db8f295]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:150151c8]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:30e40743]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:571d2d48]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:0175c56b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:180aee39]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR209C\n","truncated":false}}
%---
%[output:2e5b225f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:9557eb60]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:472a64f7]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:38bea1d6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:51193ce2]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:62e7b090]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:5eaf5562]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:3b81fabb]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:6e96a0b9]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:11c2b8d8]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:271e9aa0]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:2e868d11]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:50990ed8]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:48a070de]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:49b57f06]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: HO4P(m_Nc)  ->  HO4P\n","truncated":false}}
%---
%[output:0377cd9e]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:04c7a9a7]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:47842898]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:285983e7]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:8227a011]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:279d5232]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:9f173541]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.1\n","truncated":false}}
%---
%[output:678bb670]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:8a647cb4]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0511         | glycogen phosphorylase\n","truncated":false}}
%---
%[output:4004b77c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:47c2e24a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:4a2d96df]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:6b48ece5]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:8a3b1f01]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:2d76f749]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:7b216038]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:2791b88e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:3a46bd0a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:68ff350f]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:28292ae4]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:49d3fa21]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:0043c767]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:189f0744]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:2aa02b52]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:391354f6]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:5381724c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:514aab9d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:84ebb599]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:044b82c6]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:34f29fbc]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:03bd4efe]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:430a6195]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:9c342979]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:5fecbd07]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:968353fa]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:9d44bcb0]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:9350e72d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:23048000]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:3cdc39d2]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:64f1fd7c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:11ccd763]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:28583c2c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:0879dd9f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:102b482d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:8ee2eab6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:24dfe9ca]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: HO4P(m_Nc)  ->  HO4P\n","truncated":false}}
%---
%[output:36265e10]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:9d16ed42]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:51d276c2]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:6ac165b3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:3b72d254]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:2b43278c]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:38527866]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:79160696]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:71f17d2d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:30085f2c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:6caf7383]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:4dba2507]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:04bacea8]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR209C\n","truncated":false}}
%---
%[output:060e2f3b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:5210269f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:6278097a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:5ef7d6e3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:613cb2ac]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:73d7674c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:212e233a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.7.7.5\n","truncated":false}}
%---
%[output:8da31ad0]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:18d4a4f6]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:73a5baed]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:02807010]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:6b9a8477]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:8c5a63e6]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:85721dcd]
%   data: {"dataType":"text","outputData":{"text":"\n--- Metabolites ---\n","truncated":false}}
%---
%[output:707dc87e]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:049a7b31]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:4acad47d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:973a0b24]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:9ef2949a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:225af4b3]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:4d224d46]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:90aa1049]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:20bd383f]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:6cd288f9]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:6c3f6104]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:80381dff]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:60d9e700]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:576a001a]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:72c81ce2]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:814cfcfe]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:9aefd86f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:7b39c395]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:18158001]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H13N5O4 + HO4P(m_Nc)  ->  C5H5N5 + C5H9O8P\n","truncated":false}}
%---
%[output:33779621]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:6a3604ed]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:7bcfb366]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:11a3ea93]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:208ceabf]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:6e39a894]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:55c447bc]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:3c3ecbb7]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:32b8b06e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:2ac7d9db]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:9bb6c6ce]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:7cab1ecf]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:31364dae]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0943         | purine-nucleoside phosphorylase\n","truncated":false}}
%---
%[output:40150084]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:078bed90]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:6e621266]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:54f1bb70]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:18e69d18]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:50986a14]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:41ffd83f]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:429405d7]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:89eaac4b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:5fb97a03]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:77cd6f4c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:54482064]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:66a1b145]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:23afbc94]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:8b04fa95]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:282295bd]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:3baedd67]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:3200e3f0]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:862c6a68]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:09c5fa8b]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:6e78db42]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:11cbaebe]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:9c9de3bb]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:7b5ec38f]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:832cd459]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:71a7d17e]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:5a93cbf3]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:1265d684]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:7c73bf20]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:2be1a440]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR209C\n","truncated":false}}
%---
%[output:251232ef]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:710f1a63]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:26bb3397]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H13N5O4 + HO4P(m_Nc)  ->  C5H5N5O + C5H9O7P\n","truncated":false}}
%---
%[output:79004910]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:34854dae]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:24b77572]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:052a4dc1]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:660c30ee]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:2beb4b35]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:807e368f]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:0479c72f]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:0b8cc50e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:9583328a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:3216d79b]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:442b20b2]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:8b96967e]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:919998b1]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:30c6d3b5]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:2dcb0013]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:3b187a08]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:703c4507]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:36e7dc5c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:3fae533b]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:20840e07]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:28b8e5df]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:065f1fdd]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:996c34ce]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:5c3258a3]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:84c73fa4]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:916ea3b1]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:77c88e0d]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:356a1c5d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:15e108a3]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:99d6e468]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:3963f35b]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:7004f370]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:0871038a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:37a54a9e]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:16012c06]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:792c2b91]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:07198641]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:22066f68]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:7c55697a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:696bda48]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:80227d2a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:92595dc1]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:9039aeab]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:5840ad9d]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | NADP(+)                        | C21H25N7O17P3\n  m_d   (differentiator output \/ counter input)    | 2-oxoglutarate                 | C5H4O5\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | L-glutamate                    | C5H8NO4\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | phosphate                      | HO4P\n","truncated":false}}
%---
%[output:78747869]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + HO4P(m_Nc)  ->  H + HO4P\n","truncated":false}}
%---
%[output:7a3cc764]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:9997b74e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:660829cf]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:00d557cf]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | ammonium                       | H4N\n  m_d   (differentiator output \/ counter input)    | L-glutamine                    | C5H10N2O3\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | phosphate                      | HO4P\n  m_Cs  (counter cofactor)                         | glyceraldehyde 3-phosphate     | C3H5O6P\n  m_Nc  (NN convergence output)                    | adenine                        | C5H5N5\n","truncated":false}}
%---
%[output:5a9a9f43]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | L-glutamate                    | C5H8NO4\n  m_Df  (fast intermediate)                        | ammonium                       | H4N\n  m_d   (differentiator output \/ counter input)    | L-glutamine                    | C5H10N2O3\n  m_Ds  (slow branch product)                      | L-aspartate                    | C4H6NO4\n  m_c   (counter output \/ NN input)                | phosphate                      | HO4P\n  m_Cs  (counter cofactor)                         | glyceraldehyde 3-phosphate     | C3H5O6P\n  m_Nc  (NN convergence output)                    | 5,6-bis(diphospho)-1D-myo-inositol tetrakisphosphate | C6H7O30P8\n","truncated":false}}
%---
%[output:59d108ac]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:873a3aa9]
%   data: {"dataType":"text","outputData":{"text":"  m_in  (differentiator input)                     | 10-formyl-THF                  | C20H21N7O7\n  m_Df  (fast intermediate)                        | 5'-phosphoribosyl-N-formylglycineamide | C8H13N2O9P\n  m_d   (differentiator output \/ counter input)    | L-glutamate                    | C5H8NO4\n  m_Ds  (slow branch product)                      | THF                            | C19H21N7O6\n  m_c   (counter output \/ NN input)                | 2-oxoglutarate                 | C5H4O5\n  m_Cs  (counter cofactor)                         | ammonium                       | H4N\n  m_Nc  (NN convergence output)                    | 2-oxoglutarate                 | C5H4O5\n","truncated":false}}
%---
%[output:8ca8d41f]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:37a5c6f0]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:34bcf3b9]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:0307d386]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:64e11a91]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:758609e1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:8bbb5994]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:10dd463b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:5d38b478]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0949         | guanosine phosphorylase\n","truncated":false}}
%---
%[output:3e4638ac]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:7a1b5f28]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:850b854d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:17fe2835]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:3ead8dbf]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:4a408b36]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:422c3da0]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:48d3432a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:67010e23]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YCL050C\n","truncated":false}}
%---
%[output:11aed6ee]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:856c7768]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:6018e26b]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:8361fbbc]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:0470a6f3]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:57bb0fae]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:486272e7]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:00a9ceb4]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:5afbe187]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:7efec691]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:371021aa]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:5c16f8ed]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:721e3de0]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:5eb400bf]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:108242b7]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:8a0461ad]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:83380371]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:4906a666]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:7426e05e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:793bfa7c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:982f265d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C11H15N5O3S + HO4P(m_Nc)  ->  C5H5N5 + C6H11O7PS\n","truncated":false}}
%---
%[output:0f392906]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:69737e47]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:0597eebb]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:87bf322e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:3d7f97ac]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:80512382]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:3b28b432]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:0b74322e]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:65af014d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:3ea0802b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:5c95cfb2]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:7beb33cc]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.1\n","truncated":false}}
%---
%[output:3d22e6f9]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:929e1ca3]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:29d4a32c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:5c44038e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:5bbd6df2]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:43b88841]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:3458616f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:8bf9456d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:6daf15a0]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:1ffb7800]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:33b49a8b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:10766e25]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:63e4edbd]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:60e04d4b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:5f38dcd7]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:4bcc4f48]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:511d1ab9]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:44324e49]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:8cbcd850]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:80a42bb2]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:97509e88]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:31e83f45]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:1421e4c1]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:39ee2afe]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:9463d9f6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:477bdec0]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:80efb81f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:0003fe7a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:698fe15b]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:9b981943]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:1b71d38a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:7edff9fd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:2b59ddc7]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:78e0d8cb]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H6O27P7 + 2 H + HO4P(m_Nc)  ->  C6H7O30P8 + H2O\n","truncated":false}}
%---
%[output:3bf4a395]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:2d907fe9]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:21d1fe90]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:9d63ae74]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:76cefc60]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:97e0b1a4]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:5b4c6133]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:1cfbb3ea]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR013C\n","truncated":false}}
%---
%[output:0b61c299]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:4ca37162]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:44ffcf62]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:41eac6d3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:6be415fa]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:71ab5a36]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:0a599559]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:3402fc27]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:6f407383]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:5879e098]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:8bd0324a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:3c27f940]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:6723cec2]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0951         | inosine phosphorylase\n","truncated":false}}
%---
%[output:4d930ac9]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:546a17f8]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:8a63e90b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:4d7dd10d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.1\n","truncated":false}}
%---
%[output:45ef46e5]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:5a162b5c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:65533ffa]
%   data: {"dataType":"text","outputData":{"text":"\n--- Reactions ---\n","truncated":false}}
%---
%[output:60638f3f]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:35269413]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:94fd2bb3]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:542eb6f1]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:01f0a659]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:05d3d0cc]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:2a39eba0]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:7f7675e6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:9557989e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:2a2f66e1]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:588e3100]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:6dd39b5c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:8cf84f45]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:4870c334]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:78380565]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:9c04d3c6]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:14bd5abd]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:0a613f40]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H6O27P7 + 2 H + HO4P(m_Nc)  ->  C6H7O30P8 + H2O\n","truncated":false}}
%---
%[output:1cfd6e48]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:54cf6d47]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_1026         | sulfate adenylyltransferase (ADP)\n","truncated":false}}
%---
%[output:233ad83c]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:3ce4967e]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:65a2ad73]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:8cf99761]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:7c100058]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:652b04c4]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:109503f7]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:8d4725f1]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:7b77b78e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:55e04666]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:288558c7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:5622632e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:2c2f7a22]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:124de497]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:91d3d3b7]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:71a83d36]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:5cfe118a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:005dda0f]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:562d04da]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:5954ce2b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:324b21e9]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:59132e44]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_2008         | phosphate transport\n","truncated":false}}
%---
%[output:3212ffb6]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:257d2b8a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:1510e382]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:34b4676c]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:8bbc8691]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:9ca9b3b1]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:2785b065]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:6f44788f]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:5b44c6cb]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:017b4ce0]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:987e5db1]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:592173d2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:6ea134fc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:3c18bc7a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:92701816]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:5a93c114]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR013C\n","truncated":false}}
%---
%[output:20bb04af]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:0a731811]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:7e79642e]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:0db64fd5]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:7b2d376c]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_3649         | phosphate transport, cytoplasm-vacuolar membrane\n","truncated":false}}
%---
%[output:2b86b152]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + 2 H + C6H7O21P5 + 2 HO4P(m_Nc)  ->  C6H7O30P8 + C10H12N5O10P2 + 2 H2O\n","truncated":false}}
%---
%[output:7f54a3b5]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:5a1ee8ff]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:34f8cc25]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:7d274839]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:9911cb42]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:53cd97f9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:7ad2409a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:71598952]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:6aa21168]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:9c723b70]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:7fbba03d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:1c671c66]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:98161b69]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:8144d6c2]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:0f2eed59]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:1325f9a4]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:28ecd2a6]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_4749         | adenosine phosphorylase\n","truncated":false}}
%---
%[output:13f86f3e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:8f0fa6a8]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:1f44d949]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:73a953e1]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:9629fa64]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:92d1bcdb]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:5f01f85a]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:80ced352]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:7d25b0e0]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:65193b44]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:2a5edc41]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:1e09b7a7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:5b45c4af]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:1becb455]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:2b2cb234]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:18a9961c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:77245c68]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:83cabba2]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:9218fb29]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:5a86ac05]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:874b7acd]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:7b789579]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:9536904f]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_4758         | deoxyguanosine phosphorylase\n","truncated":false}}
%---
%[output:55b6f093]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:6a607bda]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:6b7db9ee]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:6c459ae6]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:18723b21]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:8b87cd21]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H10O5 + HO4P(m_Nc)  ->  C6H11O9P\n","truncated":false}}
%---
%[output:348604a6]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:14292986]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:299de291]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:97448a92]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:0e735326]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.28\n","truncated":false}}
%---
%[output:47700dae]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:482bbb64]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:627299b7]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:3aba06ee]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:8cb335a2]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:1264200e]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:8465fe82]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:74660c87]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:449791b6]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:941ddc56]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_1245         | phosphate transport\n","truncated":false}}
%---
%[output:8070cf7e]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:32655eef]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:3d541ac5]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:156b3164]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:4a71bad3]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR209C\n","truncated":false}}
%---
%[output:308660f7]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:693746f5]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:38a5c276]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:8453bb16]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:25950ec6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:74f69307]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:36b2fdcd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:23eb441a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:09578b1c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:3bdadd49]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:17f95153]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:756cb904]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:8aff8459]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:615f906d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:690a08c6]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0075         | 5'-methylthioadenosine phosphorylase\n","truncated":false}}
%---
%[output:7cd56229]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:134c19ad]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:853afc98]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:12af92de]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:7c6daa0f]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:888387ac]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:1b4a6c9e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:055cfacb]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:5ab04c2d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:84ca39c9]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:77f4d03a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C11H15N2O5 + HO4P(m_Nc)  ->  C5H9O8P + H + C6H6N2O\n","truncated":false}}
%---
%[output:141114fb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:64b5f075]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:6469fd3c]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:772e76d7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:08e869d2]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:60028875]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:26559a5c]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:7a3d0a7c]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:60116d6e]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0088         | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:90999d2f]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:0ca3fbb0]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:6156936f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:4d561a93]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:29e21551]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:4adef79a]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:2c74ca67]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:706c5125]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.7.4.21;2.7.4.24\n","truncated":false}}
%---
%[output:1fa89868]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:13eb50e5]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:825896b9]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:2bfe58f1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:9a3a39e5]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:96f30d5d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:63034c2d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:28327f8d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:35939cee]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:5b17f74e]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:9e82fe58]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:35e6ed0d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:595cade5]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:341926ac]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:7229a42b]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0089         | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:0e9e6037]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:1ffefa54]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:8bc3b2f3]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:31ab031b]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:38e22857]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:01bcadbe]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:2aa4393d]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:8e408e64]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:6f34d664]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:8c150137]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:529926c4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:76d1ff4c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:56ef8566]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:4d74db61]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H13N5O5 + HO4P(m_Nc)  ->  C5H9O8P + C5H5N5O\n","truncated":false}}
%---
%[output:545b86d3]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:510adde2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:2619dac1]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:0ff528b6]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR209C\n","truncated":false}}
%---
%[output:5b6a8933]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:59de3bd1]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:20e19cdc]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0358         | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:491bd3f1]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:79f9df68]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:1b67ca91]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:877bbda3]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:83bea097]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:225fa24c]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:553f6c3a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:85b2da6e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:6a8f15d3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:723fcdbf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:3403f4bb]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:29517250]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:89a00686]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:2c6692ac]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:6c65914e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:415f3196]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:7829b937]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:06e5ffcc]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:0ccdf30e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:4dae2b85]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0511         | glycogen phosphorylase\n","truncated":false}}
%---
%[output:809bd8a1]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:8f5fd96c]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:979e90c9]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:3019662b]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.7.4.21;2.7.4.24\n","truncated":false}}
%---
%[output:811af5b5]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:042758ba]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:4a0039cd]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:2b2cc5ee]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:0a1727a4]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:3c90121d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:34492b08]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:9431c610]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:314f7d37]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:18ed485d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:6b95b09f]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:1b7466a0]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:6265b461]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:045f3b73]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:9f596d47]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N4O5 + HO4P(m_Nc)  ->  C5H9O8P + C5H4N4O\n","truncated":false}}
%---
%[output:3f27e2be]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:5dd6406b]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0943         | purine-nucleoside phosphorylase\n","truncated":false}}
%---
%[output:5e6d15eb]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:09661b1f]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:63346b63]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:33e12251]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:711ed926]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:2cf2746d]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:0ad804bd]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:100ab65d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:68c00048]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:176a1a49]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:49c92ec8]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:494a66fe]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:547544ba]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:7337ec59]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:6b4064a3]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:631c4891]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:2c713c3d]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:7bcd1b99]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:1f790cba]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:3931a7fe]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:12abdea9]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:7ff0ca6e]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0949         | guanosine phosphorylase\n","truncated":false}}
%---
%[output:59d62410]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:45c6b1cb]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:895feaf5]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:6f4d3815]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:624dce73]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:3bcc992f]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:100a3259]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:6cc20763]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:8d49acd7]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:084ae19c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:2539792e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YER053C or YJR077C\n","truncated":false}}
%---
%[output:1e475c9d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:7f69e7a8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:0e7bf32f]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:670a661f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:0ac58742]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:3239fec5]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:77391a3d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:90345e89]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:5506aaa0]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.7.4.21\n","truncated":false}}
%---
%[output:1c2ff48c]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0951         | inosine phosphorylase\n","truncated":false}}
%---
%[output:5d80bbf9]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O10PS + HO4P(m_Nc)  ->  C10H12N5O10P2 + H + O4S\n","truncated":false}}
%---
%[output:81ca3f42]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:7c55fbae]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:56442f23]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:00cd6647]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:49ed3f3a]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:1e6deb54]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:17203297]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:197774a5]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:4f0de932]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:718dfdf0]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:1a467ed6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:438d7e38]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:563c053e]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:3107230c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:52c31584]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:6ee502fa]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:6776a85c]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_1026         | sulfate adenylyltransferase (ADP)\n","truncated":false}}
%---
%[output:1c7d0743]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:81154071]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:72a8caef]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:2309ebd3]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:756b2e28]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:208e3c33]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:0d57fddd]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:1dfaddbd]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:066cb8ca]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:733cf8df]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:1ee35c4a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:2cf198fe]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:639977bf]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:3015739a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:709e25ac]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:4dbb3e9d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:9c09908d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:7c200491]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:19230c7f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:3168c4a7]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:198ba20e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:4d3a21a3]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:64529e4f]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_2008         | phosphate transport\n","truncated":false}}
%---
%[output:4ebac246]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:2983c80a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:00aacdba]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:60eedb4d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:8ff21982]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:264c9950]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: HO4P(m_Nc)  ->  HO4P\n","truncated":false}}
%---
%[output:0782c66e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:2243b79f]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:805fd59c]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:9967f689]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:72dd4bab]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:54723426]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:938f5a07]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:03afeb76]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:50e4f630]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.1.1\n","truncated":false}}
%---
%[output:94d31078]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:94b5812c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:87fc5a05]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:2ea93c13]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:8ac45805]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:01e41403]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_3649         | phosphate transport, cytoplasm-vacuolar membrane\n","truncated":false}}
%---
%[output:923ba3d0]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR017W\n","truncated":false}}
%---
%[output:0972ae65]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:22459553]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:787e4ca6]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:6cfb7f07]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:158b2327]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:98247548]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:9559d325]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:82af6900]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:1ec93e3f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:52ee7670]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:0c62363d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:93d01694]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:629ab37e]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:33300347]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:28417d3e]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:8574ef95]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:77164893]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:4835e076]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:5eac50b9]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_4749         | adenosine phosphorylase\n","truncated":false}}
%---
%[output:6c12301f]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:4a2c42e0]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:3a52db53]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:27ec8cd7]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:52b6686d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:54054711]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:88d49548]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:70c99a1e]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:4fd5fb83]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:82c0d72f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:1eb13aa2]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: HO4P(m_Nc)  ->  HO4P\n","truncated":false}}
%---
%[output:247a78a0]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:23fe02c5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:2b10ba8d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:13ca59f2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:9b88ed30]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:4fda50c6]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:4b2c9e37]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:74955944]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:42ad7701]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_4758         | deoxyguanosine phosphorylase\n","truncated":false}}
%---
%[output:55512952]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:39591b81]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:5ee33073]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:56a7931e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:55e57b72]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:69cdad03]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:9b897334]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:078c14be]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:6c65d7d7]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:3140324a]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:6dd38e0e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:3bbe9e20]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:6285a20c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:4d3d7188]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.1;2.4.2.28;3.2.2.3\n","truncated":false}}
%---
%[output:5f665678]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:4c483888]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:45103c3f]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:4c66c717]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:462430cc]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:8039e85a]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:363dc789]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:342fdf63]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0075         | 5'-methylthioadenosine phosphorylase\n","truncated":false}}
%---
%[output:8a826441]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:32e401f6]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:880eba82]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:7c109e76]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:32874822]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:02118ba6]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:5b385c24]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:7f5812c7]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:835d07a6]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:420d3b5d]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:5a0b3f59]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:783dfb41]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:1ee3b958]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:9d2222dd]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR410W\n","truncated":false}}
%---
%[output:3eb6822a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H13N5O4 + HO4P(m_Nc)  ->  C5H5N5 + C5H9O8P\n","truncated":false}}
%---
%[output:067e16a2]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:8bf81094]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:5a303777]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:4b75ddbe]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:29e37ae4]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:73055ee5]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0088         | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:19948de7]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:75beb517]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:6be3e630]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:849c2831]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:70b7bcc0]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:31f81824]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:9cebd55f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:939fd7d4]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:649dac2b]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:77a02292]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:957104f1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:06a56b46]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:14305eaa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:28c1c108]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:61054bab]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:2a1a733b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:370fc088]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:0475ad07]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:5caab718]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:547c88d6]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:649dc9e9]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0089         | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:5e5e25b1]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:1effe680]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:166daa85]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:12d08640]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:958144ad]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:2b5887ca]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:29d6354f]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:1eaa5ed2]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:98c1c0cc]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.1\n","truncated":false}}
%---
%[output:404dc7ff]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:915d5299]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:17fe234c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:48c595be]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:8382b765]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:20160692]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:57374d5b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:9b22f22c]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:48799ffa]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:582ea795]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H13N5O4 + HO4P(m_Nc)  ->  C5H5N5O + C5H9O7P\n","truncated":false}}
%---
%[output:892431f2]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0358         | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:2610d67b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:6de7f332]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:0c332368]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:60f0a685]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:184d77d5]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:568a297a]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:37cec4d3]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:35fa9a7b]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:4ad23462]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:9ff3da68]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:052586da]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:456c5c69]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:42735cf0]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:3e762204]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:34007a44]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:3179190b]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:7bedcd25]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:58cccef7]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:6357edde]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:39d69461]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:894d3deb]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:0498429a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:15acaac5]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0511         | glycogen phosphorylase\n","truncated":false}}
%---
%[output:70706130]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:10286563]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:8abd0b2c]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:2eb48a50]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:1840d493]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR410W\n","truncated":false}}
%---
%[output:7542eeac]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:0ed45ef0]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:00f07ee0]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:0db968ab]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:3a22683e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:52b7a64d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:17f40f08]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:1623b9d5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:0723ccaa]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:96bdb53e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:91ce738c]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:174c1948]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:8dd7d69d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:6922ac8b]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0943         | purine-nucleoside phosphorylase\n","truncated":false}}
%---
%[output:183c372d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:0baccec1]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:3dce8aa4]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C11H15N5O3S + HO4P(m_Nc)  ->  C5H5N5 + C6H11O7PS\n","truncated":false}}
%---
%[output:997db68b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:8550f37b]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:6fd57591]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.1\n","truncated":false}}
%---
%[output:2aeda85d]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:01058bcc]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:71dbe905]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:4fb4b327]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:19903fdb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:2894550b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:00d4f58c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:4208a70d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:8514f1bb]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:6721459e]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:4988c375]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:01cacfb3]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:32663a74]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0949         | guanosine phosphorylase\n","truncated":false}}
%---
%[output:3c66d2dd]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:81540859]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:5542d5da]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:36e9588b]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:614d4d7f]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:3722b70b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:7f773905]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:637f710d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:78d880ab]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:7337673a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:3cd9242c]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:3eb168dd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:6ac5d2fd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:9982a720]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:3d0cb05a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:1f9ceee8]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:20afdda3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:93014233]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:2d28c47a]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:85e5f7b2]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:3e9a918a]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:904c8642]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:083d58d6]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:44ba9a0e]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0951         | inosine phosphorylase\n","truncated":false}}
%---
%[output:5e838c69]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:0d758d2d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:25cab39a]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:96f00851]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:6032c1a5]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:978cbff1]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H6O27P7 + 2 H + HO4P(m_Nc)  ->  C6H7O30P8 + H2O\n","truncated":false}}
%---
%[output:40cda491]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:8c5ab970]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:5bd591aa]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:8f89d508]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:0a1d04ff]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:6590dd0c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:1a285fde]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:1976301d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:31de50d8]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:25ddbc36]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:9a8397da]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:0c56f190]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YDR017C\n","truncated":false}}
%---
%[output:393c1bb0]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:15caa544]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:348d9cba]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.7.7.5\n","truncated":false}}
%---
%[output:24fbfaba]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_1026         | sulfate adenylyltransferase (ADP)\n","truncated":false}}
%---
%[output:7968c866]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:44e4dd66]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:37f867fa]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:2ae8b8e9]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:562c2813]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:9519a2c5]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:0a97bc81]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:27701269]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:06bfb699]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:8c986a45]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:5d40f319]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:67944c16]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:484a2858]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:417226df]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:2d67fd90]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:8c7ff89d]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_2008         | phosphate transport\n","truncated":false}}
%---
%[output:4c88cd3b]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:8acacfdf]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:3c271d5b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:8ba440b2]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:87e3c15b]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:726d575b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:44a36b77]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:8d7afd63]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:531dfbc1]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:724a0098]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:5dd21c1e]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:5d756a61]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:9f1d5ebe]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:087d6db6]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H6O27P7 + 2 H + HO4P(m_Nc)  ->  C6H7O30P8 + H2O\n","truncated":false}}
%---
%[output:3e07d1a3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:6bbfb28f]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:4bc6d90e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:95e20770]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:1c52fc70]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:56b44529]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:48461f64]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:9f12b557]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:18b1c781]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_3649         | phosphate transport, cytoplasm-vacuolar membrane\n","truncated":false}}
%---
%[output:98aed323]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:219952e6]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:8707e8dc]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:248be4f5]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:27209e67]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:17f12838]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:64b95e3d]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:7dec6592]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:9912d86a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:7a3a74e7]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:7c8fb0c0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:4a6017fe]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:2be5abb0]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:809c06c9]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:4a1d01d9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:31ededcf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:9b49b12a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:99a1bddd]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:0d69a77e]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:7223afa0]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:7e59b6c1]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:86d61545]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:84603a2c]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_4749         | adenosine phosphorylase\n","truncated":false}}
%---
%[output:10278aa8]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:122cc42c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:2ba89321]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:41aedba7]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:6d075e35]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:6abad913]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:11b755ab]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:38662d67]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR160W\n","truncated":false}}
%---
%[output:62e51676]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:9cb00fba]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:56175571]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:67b68563]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:49a848b3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:2780cd19]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + 2 H + C6H7O21P5 + 2 HO4P(m_Nc)  ->  C6H7O30P8 + C10H12N5O10P2 + 2 H2O\n","truncated":false}}
%---
%[output:1e79beea]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:1ea8f1f9]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:83a146cd]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:1e149a03]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_4758         | deoxyguanosine phosphorylase\n","truncated":false}}
%---
%[output:68dcd6a2]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:9f929bf7]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:500d110d]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:6e2e6f98]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:0e139b41]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:9afbdb24]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:25f02230]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:6f44d599]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:1efffd2a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:9328dd14]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:972cb268]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:13b16224]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:17354a54]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:202c22cf]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:235b6b80]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:7dee650a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:76e886ba]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:10ca7a85]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:2874229c]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:98236aac]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:53265cb6]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:6f6bd5ec]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:281445fb]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_3605         | phosphate transport, cytoplasm-cell envelope\n","truncated":false}}
%---
%[output:07818ad6]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:4b8db7db]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:9d260a07]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:2f836a3c]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:2afff1a9]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:77db4229]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:7f972642]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:2c4f9426]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:1279f925]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:13535b23]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:5fbfd73e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:65a42c56]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:1092e82c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:1e79abd1]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:2a223d6d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:6546d9dc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:936e8057]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:55200965]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:4f97d10b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H10O5 + HO4P(m_Nc)  ->  C6H11O9P\n","truncated":false}}
%---
%[output:1c9db429]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:280dc209]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_1245         | phosphate transport\n","truncated":false}}
%---
%[output:37a41491]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:03cfe47a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:450b4592]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:6cbc338b]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:0871c531]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:2ff5844b]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:721cc78e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:6b733f1f]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:758334d7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:711de2f7]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:0aaf3c01]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:6c2758f6]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:714c5380]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:422ab7f2]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:849fd9f2]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:13093f51]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:98c39a44]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:27647209]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:7b6eb432]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_1245         | phosphate transport\n","truncated":false}}
%---
%[output:87df82db]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:25187322]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:8487a50e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR209C\n","truncated":false}}
%---
%[output:16381bdf]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:15df4db2]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:6dc67b30]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:7a2db19e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:332bb081]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:3d10e408]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:63a01f3d]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:0696a6a5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:1f449b17]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:3fce5783]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:96de334a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:501d9f5f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:4db06b73]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0346         | dihydrofolate synthase\n","truncated":false}}
%---
%[output:83b85de5]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:470f987c]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:57009c43]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:2d11f1cf]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:3f1778a2]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:20581251]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_3605         | phosphate transport, cytoplasm-cell envelope\n","truncated":false}}
%---
%[output:2dbf6996]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:6a07e951]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0989         | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:3df9d41e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C11H15N2O5 + HO4P(m_Nc)  ->  C5H9O8P + H + C6H6N2O\n","truncated":false}}
%---
%[output:905f6834]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:4d949db6]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:8f5924d2]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:9202b44f]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0659         | isocitrate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:2d7f544d]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0216         | aspartate transaminase\n","truncated":false}}
%---
%[output:9f59faa4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0068         | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:5fe95c90]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.1\n","truncated":false}}
%---
%[output:987a3b13]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:785dd518]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:28242d5a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:74945751]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:094181d6]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0460         | gamma-glutamylcysteine synthetase\n","truncated":false}}
%---
%[output:3d762783]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:2eced40b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:220fc470]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:6706c0de]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_3605         | phosphate transport, cytoplasm-cell envelope\n","truncated":false}}
%---
%[output:2dc8503b]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0470         | glutamate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:998cfb83]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:0de969d0]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:61bc45ec]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:76f6306a]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:2df77f7d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:3f057189]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:82da4b83]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0211         | asparagine synthase (glutamine-hydrolysing)\n","truncated":false}}
%---
%[output:0200f8c6]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:621cb511]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:27bb10c3]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:211e9d7a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0079         | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:25ab0747]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_4211         | D-ribose 5-phosphate,D-glyceraldehyde 3-phosphate pyridoxal 5-phosphate-lyase (glutamine-hydrolyzing)\n","truncated":false}}
%---
%[output:7eec5c98]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:4b9efe65]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:7b1acf87]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0486         | glyceraldehyde-3-phosphate dehydrogenase\n","truncated":false}}
%---
%[output:469a23df]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:5a754beb]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:21452143]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0075         | 5'-methylthioadenosine phosphorylase\n","truncated":false}}
%---
%[output:6fe38be9]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:30fe8ec8]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_4749         | adenosine phosphorylase\n","truncated":false}}
%---
%[output:6b15f1bf]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:9af411f2]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0138         | adenine deaminase\n","truncated":false}}
%---
%[output:0fd2714b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:671a4fe6]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:188bb068]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:2cecce94]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0470         | glutamate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:0878a074]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_1118         | aspartate-glutamate transporter\n","truncated":false}}
%---
%[output:7566d7a7]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:4697bd6a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H13N5O5 + HO4P(m_Nc)  ->  C5H9O8P + C5H5N5O\n","truncated":false}}
%---
%[output:8490d017]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:72920c73]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:2d99881a]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_0211         | asparagine synthase (glutamine-hydrolysing)\n","truncated":false}}
%---
%[output:2b78ba07]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0079         | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:453716cd]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:02257bf8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_4211         | D-ribose 5-phosphate,D-glyceraldehyde 3-phosphate pyridoxal 5-phosphate-lyase (glutamine-hydrolyzing)\n","truncated":false}}
%---
%[output:279ee2a2]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR209C\n","truncated":false}}
%---
%[output:04597955]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0486         | glyceraldehyde-3-phosphate dehydrogenase\n","truncated":false}}
%---
%[output:01377b30]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:3a27bc2d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:6c5c0671]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_0089         | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:268bbce2]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_0358         | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:1c2e2eb3]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0357         | diphosphoinositol-1,3,4,6-tetrakisphosphate diphosphohydrolase\n","truncated":false}}
%---
%[output:32cba3d9]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:057efccc]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:3834a53c]
%   data: {"dataType":"text","outputData":{"text":"  r_DF (differentiator: fast)      : r_0499         | glycinamide ribotide transformylase\n","truncated":false}}
%---
%[output:95e395df]
%   data: {"dataType":"text","outputData":{"text":"  r_DS (differentiator: slow)      : r_0695         | L-tyrosine N-formyltransferase\n","truncated":false}}
%---
%[output:92a08e0d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:1b349ca6]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp (differentiator: producer)  : r_0079         | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:9d7b33a5]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.1\n","truncated":false}}
%---
%[output:9e3489df]
%   data: {"dataType":"text","outputData":{"text":"  r_Di (differentiator: integrator) : r_1031         | tetrahydrofolate:L-glutamate gamma-ligase (ADP-forming)\n","truncated":false}}
%---
%[output:0dead9ec]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:9444c103]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp (counter: producer)         : r_0018         | 2-aminoadipate transaminase\n","truncated":false}}
%---
%[output:07119c52]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:53086a1e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg (counter: gate)             : r_0476         | glutamine synthetase\n","truncated":false}}
%---
%[output:00827daf]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:169dd99a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr (counter: reset)            : r_0471         | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:1a677094]
%   data: {"dataType":"text","outputData":{"text":"  r_Na (NN: neuron a)              : r_1099         | 2-oxoadipate and 2-oxoglutarate transport\n","truncated":false}}
%---
%[output:8f5711ec]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:1681755e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:4e25ffd2]
%   data: {"dataType":"text","outputData":{"text":"  r_Nb (NN: neuron b)              : r_1112         | AKG transporter, mitochonrial\n","truncated":false}}
%---
%[output:4f0e8a9e]
%   data: {"dataType":"text","outputData":{"text":"  r_No (NN: merge\/output)          : r_0674         | L-alanine transaminase\n","truncated":false}}
%---
%[output:0a7e184a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:1d5a3c80]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:8f383d3d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:94d778d0]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:9d1dea25]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:92753874]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:5fe0a26d]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:65f1cf8a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:69014988]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:5f7cf004]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:53423dd2]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:2288493f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:4f2224c7]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:2c68e2fe]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N4O5 + HO4P(m_Nc)  ->  C5H9O8P + C5H4N4O\n","truncated":false}}
%---
%[output:5baecd8d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR209C\n","truncated":false}}
%---
%[output:5b29282f]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:84cd1fcb]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:2587dbeb]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:95469740]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:601aad3e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:107cc1a2]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:8cc0e72f]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:0cb49b29]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:63d3eca3]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:0a331101]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:24c5a872]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:1424d20b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:0724f1fd]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:0b676ae2]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:01369c25]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YCL050C\n","truncated":false}}
%---
%[output:86abaebb]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:66cadb2b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:67897898]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:0dc52fb3]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:5e905780]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:5bb5f5dd]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:88c38530]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:88c74c05]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:1378f98b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:8d4004b5]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:5b4fdcb8]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:7eb7849b]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:8bf8af8e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:1b27e458]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:2e00c8d8]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:0dfadf2c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:335c9cd2]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.28\n","truncated":false}}
%---
%[output:37a19b55]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:1d939c0b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR013C\n","truncated":false}}
%---
%[output:085aa13e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:96e3cfde]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:51e41887]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:2447c81b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:8f39e3ca]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:2bf83424]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:2bb54fae]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:215390cf]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:90c1500a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:787f59b8]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:1979521d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O10PS + HO4P(m_Nc)  ->  C10H12N5O10P2 + H + O4S\n","truncated":false}}
%---
%[output:954b8e95]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:79f2b9ee]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:88379c21]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR013C\n","truncated":false}}
%---
%[output:7eaf6d82]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:1559c752]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:8bde2a3a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:3302ef6f]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:9774d8ff]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:2c754931]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:55ab21fc]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:909e7f97]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:51bcb6bd]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:02652d33]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:8ee465a8]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:0089a340]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:16710877]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:584b9ff7]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:4db6c03f]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:97040cbf]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR209C\n","truncated":false}}
%---
%[output:8cf98f9f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:71fb6e2a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:9d5ac708]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:35032c67]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:4060240e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:9d88e280]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:31ff35db]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:6d29648c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:82b4c3a6]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:61ebbaf8]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:34fd466c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:310fb030]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:37533561]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:0cba1e27]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:325a5cea]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:4634113f]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:84c9d427]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:80172205]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR209C\n","truncated":false}}
%---
%[output:4442ac63]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:993cb599]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:4a2f6fb7]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:2c7ccea4]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:2398d181]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:3c44c574]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:6498b959]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:1a1a1bba]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:3b25e4ed]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: HO4P(m_Nc)  ->  HO4P\n","truncated":false}}
%---
%[output:4ee8d51e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.7.4.21;2.7.4.24\n","truncated":false}}
%---
%[output:11967c4a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:5e2f29d3]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:8a02a723]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:3fa6693a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:5d7ec33c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR017W\n","truncated":false}}
%---
%[output:9f52ed87]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:0f48f0c7]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:5490f5c1]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:353cefec]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:1d637d4c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:6701619c]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:38160bd8]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:1cacf52c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:9dd636a9]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:5825699b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:19389fa5]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:7f70873e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:286055e9]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:69c52bbe]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:16dcb73e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR410W\n","truncated":false}}
%---
%[output:8365cfde]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:7b020fc1]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:84c85ab3]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:03cb4bab]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:4d97eb59]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:9148aff4]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:5cf875aa]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:3eb600a9]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:204c64f6]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:6ac931ad]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:16e93374]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:2c0071ac]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:7fff8512]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:2a1d0815]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:9bfafbc5]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:04fa4a53]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:5663a411]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:0e117794]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:00e4d0fa]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR410W\n","truncated":false}}
%---
%[output:81fa684b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:2cb0aad9]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:21c76914]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:1a108d4a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:69b2020d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:5879c3d1]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: HO4P(m_Nc)  ->  HO4P\n","truncated":false}}
%---
%[output:2226baf7]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:083fd791]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:4941d0f0]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:3b5caaa8]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:96f09537]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:54dcc1a3]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:27332270]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:9d757aa5]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:6283f531]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YDR017C\n","truncated":false}}
%---
%[output:9cb315a1]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:53f01bed]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:83a38437]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.7.4.21;2.7.4.24\n","truncated":false}}
%---
%[output:47462cd9]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:4affa2b7]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:0038c101]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:0cf01430]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:2631ccf2]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:40639224]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:597c42aa]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:12888b9a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:61c7c548]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:4fd9c140]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:2454b840]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:0bcf8625]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR160W\n","truncated":false}}
%---
%[output:87ba3e0e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:106c3f5f]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:9526e22a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:50af5348]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:55dc1250]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:5fb140b2]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:9354a170]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:947d5dc3]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:7dcf83a7]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:5758c6dc]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:88b2f203]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:8935b9a8]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:806b348d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:4f501cef]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:1bcd02d3]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:9ac5bdef]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:14827de9]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:7cf86e65]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:44c3cd9f]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR209C\n","truncated":false}}
%---
%[output:5e255b2d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:0e723180]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:0a9bd84c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H13N5O4 + HO4P(m_Nc)  ->  C5H5N5 + C5H9O8P\n","truncated":false}}
%---
%[output:4522cd14]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:5f29262c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:09af1cbc]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:03b0f43c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:5e47c4fd]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:7bae529e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:3fc10494]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:692c8754]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:1841baca]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:523828f5]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:3c6d5b6e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:7558509c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR209C\n","truncated":false}}
%---
%[output:2960e4eb]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:98df5192]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:56f0f189]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:6f8480fc]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:0fe96742]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:39b2c8bc]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:73a3cc9a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:7d41a89c]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:3d47bc10]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:17a677a2]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:75ff1a96]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:7c9826f3]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:4f11a50f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.7.4.21\n","truncated":false}}
%---
%[output:5b36837b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:08b131b9]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:841da2ba]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:6bd51730]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR209C\n","truncated":false}}
%---
%[output:43695f90]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:86662e50]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:906e0875]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:6100f547]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:59eac599]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:32825530]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:073d68f9]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:8da1fb49]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:8f612036]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:6118c065]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:7b1e2d9e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:93645534]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:0634f861]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:69f7d401]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:1edf1d9d]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:8e9de97d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YCL050C\n","truncated":false}}
%---
%[output:1c596e5a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H13N5O4 + HO4P(m_Nc)  ->  C5H5N5O + C5H9O7P\n","truncated":false}}
%---
%[output:44f7c261]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:5404073d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:92b757ba]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:3b70c5ef]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:0cd6c570]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:5e0afd21]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:007152c3]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:109ba664]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:33a516c6]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:16615950]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:40d6cef3]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:5130378c]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:24f00a60]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:5e5b0577]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:30ae46de]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR013C\n","truncated":false}}
%---
%[output:172727bb]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:4460b27f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:45a23646]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:3a95dadc]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:63c778f6]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:0af3847f]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:576ced0c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:6e71e5a6]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:6a26ccd9]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:6466dfef]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:6ea69738]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:2b09d326]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:9a8512ca]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:0bc90a0a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:4b4f8b16]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:05c32c83]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:2c6f6cb3]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR013C\n","truncated":false}}
%---
%[output:72b86ac9]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:92c70104]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:51648625]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:698e57eb]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:010e08a9]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:103f13a7]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:5f2df76d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.1.1\n","truncated":false}}
%---
%[output:0ad12e5a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:6d2cf2f7]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:7a288143]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:1ddc4e2b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:18652103]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:4ee824fc]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:92caddf3]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: HO4P(m_Nc)  ->  HO4P\n","truncated":false}}
%---
%[output:88df9472]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:4390b76d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR209C\n","truncated":false}}
%---
%[output:7618cd5e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:1443a754]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:3cb32f7d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:4d4b6fa7]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:9ad3b3c5]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:332678dc]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:5619ac77]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:31cf8627]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:31f28b69]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:020e17cf]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:67c95a80]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:92895d9e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:76051958]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:277744ba]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:217a8be8]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR209C\n","truncated":false}}
%---
%[output:802f09e7]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:672fd99e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:6852bfda]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:4c804c69]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:7cdf7b9c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:157565a4]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:5818c3f1]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:15ae0e44]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:371a5c69]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:0ebc8bf8]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:1f706c8a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:5ed1f4b6]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:0741df89]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:20382e81]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:2a8c1dc0]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:32170d8d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:861d6497]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:50458c07]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YCR037C or YJL198W or YML123C\n","truncated":false}}
%---
%[output:40e51824]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:6f7ebc90]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:2c69ed35]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:16df5828]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:2715ee4a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:225bfc75]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:6bb095b0]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:265ee563]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:60e3950f]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:107ac611]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:68af04fa]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + HO4P(m_Nc)  ->  H + HO4P\n","truncated":false}}
%---
%[output:87eda5ff]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:3bdb9577]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:249bfa72]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:8dc9886b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:785cacc8]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.1;2.4.2.28;3.2.2.3\n","truncated":false}}
%---
%[output:09858bf1]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YER053C or YJR077C\n","truncated":false}}
%---
%[output:1651321d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:98258f39]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:188a2fda]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:1bd58581]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:01d6e0c0]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:63cf6405]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:33de7d25]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:23a9c9b5]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:7963f6ac]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:33a2e385]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:9c6c4a37]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:88a3d7fd]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:406096cb]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:31dfa6a5]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YER053C or YJR077C\n","truncated":false}}
%---
%[output:7258a7e4]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:1d084f55]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:0f2f217a]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:9347860b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:305ffae5]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:07b88109]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:78914f6b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:94c8c355]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:77f7a612]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:28427410]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:068cedf5]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:79dec1f6]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:789602f7]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:665a34aa]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:6c9085a5]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR113W\n","truncated":false}}
%---
%[output:397abd96]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:6e2c958f]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YCR037C or YJL198W or YML123C\n","truncated":false}}
%---
%[output:106a5974]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:53c66c76]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:60b994b2]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNR050C\n","truncated":false}}
%---
%[output:01526602]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:1d02c128]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:9fc7c8a4]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR174W\n","truncated":false}}
%---
%[output:4590b788]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:64452154]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR027C\n","truncated":false}}
%---
%[output:89b1df8d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + HO4P(m_Nc)  ->  H + HO4P\n","truncated":false}}
%---
%[output:6bd6dff5]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR019W\n","truncated":false}}
%---
%[output:3f4cb34a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:07398029]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:002295a1]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:54ebc1ff]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YJL101C\n","truncated":false}}
%---
%[output:4b643aee]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:0dab80da]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:5b0c963b]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YCR037C or YJL198W or YML123C\n","truncated":false}}
%---
%[output:5ccb59cd]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YDL215C\n","truncated":false}}
%---
%[output:46dc8241]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:32acb1a7]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:9c748129]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:221e87da]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:8a908e03]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:81e14dc1]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR124W or YPR145W\n","truncated":false}}
%---
%[output:9899458c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR061C\n","truncated":false}}
%---
%[output:4bbc1d6a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.1\n","truncated":false}}
%---
%[output:50d063fb]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:088e37b8]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YFL060C or YNL334C\n","truncated":false}}
%---
%[output:989075f4]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR192C or YJL052W or YJR009C\n","truncated":false}}
%---
%[output:863f713d]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR017W\n","truncated":false}}
%---
%[output:146eabd3]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:1020c36c]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR209C\n","truncated":false}}
%---
%[output:46072f23]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YNL141W\n","truncated":false}}
%---
%[output:4a72b0f8]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YDL215C\n","truncated":false}}
%---
%[output:992cfda3]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:86030b55]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR021C\n","truncated":false}}
%---
%[output:6a414d77]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:1f634206]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR124W or YPR145W\n","truncated":false}}
%---
%[output:718d5e32]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:260e5fa3]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:741f71aa]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR061C\n","truncated":false}}
%---
%[output:4c9d1731]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YFL060C or YNL334C\n","truncated":false}}
%---
%[output:31e5a39c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:7032cdf1]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR192C or YJL052W or YJR009C\n","truncated":false}}
%---
%[output:15ebee78]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:76bd3147]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR410W\n","truncated":false}}
%---
%[output:573c252c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C14H13N6O3 + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + C19H19N7O6 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:747c79cb]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YDR017C\n","truncated":false}}
%---
%[output:7aea2c68]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YOR163W\n","truncated":false}}
%---
%[output:1367f547]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:168f2244]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:65f4b7c4]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YDR408C\n","truncated":false}}
%---
%[output:624bd7bd]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YDR403W\n","truncated":false}}
%---
%[output:405485c3]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YGR061C\n","truncated":false}}
%---
%[output:72871d37]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: HO4P(m_Nc)  ->  HO4P\n","truncated":false}}
%---
%[output:89101d1c]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:54b14d93]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YKL132C or YOR241W\n","truncated":false}}
%---
%[output:5daf6930]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YER152C or YGL202W or YJL060W\n","truncated":false}}
%---
%[output:6ed27562]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YPR035W\n","truncated":false}}
%---
%[output:47629766]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:27859618]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YAL062W or YOR375C\n","truncated":false}}
%---
%[output:1fb0cd1e]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YOR222W or YPL134C\n","truncated":false}}
%---
%[output:1734c8e7]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YMR241W\n","truncated":false}}
%---
%[output:435bb8c3]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:10e7ee7d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H + C6H11NO3 + C5H8NO4(m_c) + C21H26N7O17P3  ->  H2O + C11H19N2O6 + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:3e237f4a]
%   data: {"dataType":"text","outputData":{"text":"      GPR: YLR089C\n","truncated":false}}
%---
%[output:59c5305d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:9fdb9ee0]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:60f3dea8]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_c)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:1d77369a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:817ef642]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.1\n","truncated":false}}
%---
%[output:8c3bcdb0]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:5db57258]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H5O7 + C21H25N7O17P3(m_Df)  ->  C5H4O5(m_d) + CO2 + C21H26N7O17P3\n","truncated":false}}
%---
%[output:9d91c16d]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:75b224ef]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:95fe3e0a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:9ea0c718]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H6NO4(m_Ds)  ->  C5H8NO4(m_c) + C4H2O5\n","truncated":false}}
%---
%[output:229a197c]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:094e1dfb]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:93c897b1]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:066cc06c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + C4H9NO2  ->  C5H8NO4(m_c) + C4H5O3\n","truncated":false}}
%---
%[output:50454003]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:85664dac]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:9f317986]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.7.7.5\n","truncated":false}}
%---
%[output:1755fa69]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_d) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_c) + C21H25N7O17P3(m_Df)\n","truncated":false}}
%---
%[output:4c717dce]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:71c0fefa]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:596ff493]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:24505d46]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:35dcb14d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:9f07c8e4]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:563c42b7]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:2455ce65]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C3H7NO2S + C5H8NO4(m_c)  ->  C10H12N5O10P2 + H + C8H13N2O5S + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:34ecf5de]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:8cfa271f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:6e1f2b22]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:68c6a51a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_c) + C19H21N7O6  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P(m_Nc)\n","truncated":false}}
%---
%[output:77177098]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:1cc2ec99]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:09ca66c4]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:37a6723b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: HO4P(m_Nc)  ->  HO4P\n","truncated":false}}
%---
%[output:9e47b197]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:5d1edf3e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:0c2a209f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:84093a6f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:71baf1eb]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:7cc41676]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:2e190aa3]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:31b8fb47]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:56946194]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:0571f99e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H2O + C5H8NO4(m_in) + C21H26N7O14P2  ->  C5H4O5 + H4N(m_Df) + H + C21H27N7O14P2\n","truncated":false}}
%---
%[output:66704650]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:6fe207b8]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:3412f8d0]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:2edab454]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_in)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:10d3d003]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:176ce60f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:100c4b9d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:986582c7]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Df) + C10H12N5O13P3 + C5H8NO4(m_in)  ->  C10H12N5O10P2 + H + C5H10N2O3(m_d) + HO4P(m_c)\n","truncated":false}}
%---
%[output:69f0d83d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:02f02f99]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:65a904d7]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.1\n","truncated":false}}
%---
%[output:2f6836d5]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + H2O + C4H6NO4(m_Ds) + C5H10N2O3(m_d)  ->  C10H12N5O7P + HO7P2 + H + C4H8N2O3 + C5H8NO4(m_in)\n","truncated":false}}
%---
%[output:595487bb]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:33d77f95]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:0f233cb7]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:50012c9a]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C8H13N2O9P + C10H12N5O13P3 + H2O + C5H10N2O3(m_d)  ->  C8H15N3O8P + C10H12N5O10P2 + H + C5H8NO4(m_in) + HO4P(m_c)\n","truncated":false}}
%---
%[output:72b81796]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:1efed5b3]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:4afba3ed]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:55d7c6ed]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C3H5O6P(m_Cs) + C5H10N2O3(m_d) + C5H9O8P  ->  H + 3 H2O + C5H8NO4(m_in) + HO4P(m_c) + C8H8NO6P\n","truncated":false}}
%---
%[output:3e1ea0ef]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:3f81ca1d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:0de7c66e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:2031a5ef]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C3H5O6P(m_Cs) + C21H26N7O14P2 + HO4P(m_c)  ->  C3H4O10P2 + H + C21H27N7O14P2\n","truncated":false}}
%---
%[output:94ea07b3]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.1\n","truncated":false}}
%---
%[output:1263c7ee]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:3fb0e022]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:9590f45d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C11H15N5O3S + HO4P(m_c)  ->  C5H5N5(m_Nc) + C6H11O7PS\n","truncated":false}}
%---
%[output:07548101]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:2bd6660b]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:9631a3b6]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:4a599c28]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H13N5O4 + HO4P(m_c)  ->  C5H5N5(m_Nc) + C5H9O8P\n","truncated":false}}
%---
%[output:0fe32050]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:7e6dc31d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:57fdea43]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:37951ece]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H5N5(m_Nc) + H + H2O  ->  H4N(m_Df) + C5H4N4O\n","truncated":false}}
%---
%[output:57eb0124]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:5c1c7eed]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:56da561a]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:49f6aaff]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:2d94ab1b]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:16ed10c7]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:647f8331]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:31f4acc8]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:71467855]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:64350e43]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H2O + C5H8NO4(m_in) + C21H26N7O14P2  ->  C5H4O5 + H4N(m_Df) + H + C21H27N7O14P2\n","truncated":false}}
%---
%[output:78af054d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:2e39546d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:15c6d799]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:8652767b]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C4H6NO4 + C5H8NO4(m_in)  ->  C4H6NO4(m_Ds) + C5H8NO4\n","truncated":false}}
%---
%[output:988c108e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:80ad93c2]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:6541fe03]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:4dfc6a23]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Df) + C10H12N5O13P3 + C5H8NO4(m_in)  ->  C10H12N5O10P2 + H + C5H10N2O3(m_d) + HO4P(m_c)\n","truncated":false}}
%---
%[output:5800b3bd]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:3d35514e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:630c9768]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:7117daa6]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + H2O + C4H6NO4(m_Ds) + C5H10N2O3(m_d)  ->  C10H12N5O7P + HO7P2 + H + C4H8N2O3 + C5H8NO4(m_in)\n","truncated":false}}
%---
%[output:6696a670]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:2a89fa5d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:2d7bae09]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:13fe025f]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C8H13N2O9P + C10H12N5O13P3 + H2O + C5H10N2O3(m_d)  ->  C8H15N3O8P + C10H12N5O10P2 + H + C5H8NO4(m_in) + HO4P(m_c)\n","truncated":false}}
%---
%[output:4c0c01ba]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:8b7a9749]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:82425e09]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:545fbb0d]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C3H5O6P(m_Cs) + C5H10N2O3(m_d) + C5H9O8P  ->  H + 3 H2O + C5H8NO4(m_in) + HO4P(m_c) + C8H8NO6P\n","truncated":false}}
%---
%[output:8c1ebc43]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:6e92b916]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:458c0924]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:76b6bbc7]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C3H5O6P(m_Cs) + C21H26N7O14P2 + HO4P(m_c)  ->  C3H4O10P2 + H + C21H27N7O14P2\n","truncated":false}}
%---
%[output:35873879]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:8ecab688]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:1c815e97]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.12\n","truncated":false}}
%---
%[output:8d72c36c]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H6O27P7 + 2 H + HO4P(m_c)  ->  C6H7O30P8(m_Nc) + H2O\n","truncated":false}}
%---
%[output:2e48fc2b]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:0f546a4a]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:238a16a4]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.5.1.10\n","truncated":false}}
%---
%[output:8e6eb953]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + 2 H + C6H7O21P5 + 2 HO4P(m_c)  ->  C6H7O30P8(m_Nc) + C10H12N5O10P2 + 2 H2O\n","truncated":false}}
%---
%[output:74055424]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:269d765f]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.1.1.42\n","truncated":false}}
%---
%[output:66476c8e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.1\n","truncated":false}}
%---
%[output:56cdc211]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H7O30P8(m_Nc) + 3 H2O  ->  3 H + C6H7O21P5 + 3 HO4P(m_c)\n","truncated":false}}
%---
%[output:8c09256b]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.19\n","truncated":false}}
%---
%[output:4965f7d7]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:7e4b556b]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:1993e5ce]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.2\n","truncated":false}}
%---
%[output:0b9215e5]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:4851e2bf]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:333af62d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.2\n","truncated":false}}
%---
%[output:74785caa]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:6ae14656]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:427690b2]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C20H21N7O7(m_in) + C7H14N2O8P  ->  C8H13N2O9P(m_Df) + H + C19H21N7O6(m_Ds)\n","truncated":false}}
%---
%[output:03afbbac]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.5.4\n","truncated":false}}
%---
%[output:1f1522b0]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.5.3\n","truncated":false}}
%---
%[output:8c6171b5]
%   data: {"dataType":"text","outputData":{"text":"      EC : 4.3.3.6;3.5.1.2\n","truncated":false}}
%---
%[output:615610f2]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C20H21N7O7(m_in) + C9H11NO3  ->  H + C10H10NO4 + C19H21N7O6(m_Ds)\n","truncated":false}}
%---
%[output:9a426651]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.2.1.12\n","truncated":false}}
%---
%[output:191f5dfb]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.28\n","truncated":false}}
%---
%[output:37066068]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.4.2.1\n","truncated":false}}
%---
%[output:49a2e1d1]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C8H13N2O9P(m_Df) + C10H12N5O13P3 + H2O + C5H10N2O3  ->  C8H15N3O8P + C10H12N5O10P2 + H + C5H8NO4(m_d) + HO4P\n","truncated":false}}
%---
%[output:3a5b5539]
%   data: {"dataType":"text","outputData":{"text":"      EC : 3.5.4.2\n","truncated":false}}
%---
%[output:4dd83f06]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.2\n","truncated":false}}
%---
%[output:24eecd03]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:0c4b77ac]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C10H12N5O13P3 + C5H8NO4(m_d) + C19H21N7O6(m_Ds)  ->  C24H27N8O9 + C10H12N5O10P2 + H + HO4P\n","truncated":false}}
%---
%[output:676b2c9e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:3d659633]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.5.4\n","truncated":false}}
%---
%[output:6b656174]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.5.3\n","truncated":false}}
%---
%[output:2cebc2b6]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H6O5 + C5H8NO4(m_d)  ->  C5H4O5(m_c) + C6H10NO4\n","truncated":false}}
%---
%[output:6b6ae280]
%   data: {"dataType":"text","outputData":{"text":"      EC : 4.3.3.6;3.5.1.2\n","truncated":false}}
%---
%[output:5a3ba033]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.2.1.12\n","truncated":false}}
%---
%[output:2fe86fda]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.7.4.21;2.7.4.24\n","truncated":false}}
%---
%[output:4d444e2e]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: H4N(m_Cs) + C10H12N5O13P3 + C5H8NO4(m_d)  ->  C10H12N5O10P2 + H + C5H10N2O3 + HO4P\n","truncated":false}}
%---
%[output:90c77023]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.7.4.21\n","truncated":false}}
%---
%[output:2662aed1]
%   data: {"dataType":"text","outputData":{"text":"      EC : 3.6.1.52;3.6.1.60\n","truncated":false}}
%---
%[output:0b1b1a6e]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.1.2.2\n","truncated":false}}
%---
%[output:36f25790]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_c) + H4N(m_Cs) + H + C21H26N7O17P3  ->  H2O + C5H8NO4(m_d) + C21H25N7O17P3\n","truncated":false}}
%---
%[output:292235d2]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:1b060e02]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.5.3\n","truncated":false}}
%---
%[output:4b08f57d]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.2.17\n","truncated":false}}
%---
%[output:62aaa7f7]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C6H6O5 + C5H4O5(m_c)  ->  C6H6O5 + C5H4O5(m_Nc)\n","truncated":false}}
%---
%[output:026d3b12]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.39;2.6.1.57;2.6.1.7\n","truncated":false}}
%---
%[output:9b375755]
%   data: {"dataType":"text","outputData":{"text":"      EC : 6.3.1.2\n","truncated":false}}
%---
%[output:604167cb]
%   data: {"dataType":"text","outputData":{"text":"      EC : 1.4.1.4\n","truncated":false}}
%---
%[output:712c7708]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_c) + C6H5O7  ->  C5H4O5(m_Nc) + C6H5O7\n","truncated":false}}
%---
%[output:43cb3863]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:6c799ba7]
%   data: {"dataType":"text","outputData":{"text":"      EC : (none annotated)\n","truncated":false}}
%---
%[output:31c088b3]
%   data: {"dataType":"text","outputData":{"text":"      EC : 2.6.1.2\n","truncated":false}}
%---
%[output:223239dd]
%   data: {"dataType":"text","outputData":{"text":"      Eqn: C5H4O5(m_Nc) + C3H7NO2  ->  C5H8NO4 + C3H3O3\n","truncated":false}}
%---
%[output:6288ed4c]
%   data: {"dataType":"text","outputData":{"text":"\nSaved nn_full_motifs.mat (48 full motifs, ranked by gene complexity)\n","truncated":false}}
%---
