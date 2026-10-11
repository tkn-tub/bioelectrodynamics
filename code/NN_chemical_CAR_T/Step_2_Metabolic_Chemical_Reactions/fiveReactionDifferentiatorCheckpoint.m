function checkpoint = fiveReactionDifferentiatorCheckpoint(model,motifs,diagnostics,paths,settings,info)
% Store numeric motif columns once. Descriptions/equations are dictionaries
% indexed by model index, avoiding repeated text on 100000 result rows.
checkpoint = struct;
checkpoint.TopologyVersion = "five_reaction_differentiators_v2";
checkpoint.Stage = settings.stage;
checkpoint.KineticsValidated = false;
checkpoint.TopologyDescription = "Block 1: X->M1->Y1, X->F1, F1+Y1 removal, F1 removal; Block 2: Y1->M2->Y2, Y1->Y2, M2+Y2 removal, Y2 removal.";
checkpoint.OutputDescription = "m_d1 = m_d = Y1; m_F1_delay = F1; m_d2 = Y2 (sole counter signal).";
checkpoint.ReactionRoles_Block1 = ["r_S1","r_S2","r_F","r_I","r_R"];
checkpoint.ReactionRoles_Block2 = ["r_S1_2","r_S2_2","r_F2","r_I2","r_R2"];
checkpoint.CheckpointStorageVersion = 2;
checkpoint.EC_AnnotationsRequired = true;
checkpoint.EC_ModelMetIDs = string(model.mets(:));
checkpoint.EC_ModelRxnIDs = string(model.rxns(:));
checkpoint.InputAdapterPaths = paths;
checkpoint.SearchSettings = settings;
checkpoint.SearchInfo = info;
[~,codes] = getReactionECAnnotations(model);
names = motifs.Properties.VariableNames;
metRoles = names(startsWith(names,'m_'));
rxnRoles = names(startsWith(names,'r_'));
m = unique(reshape(table2array(motifs(:,metRoles)),[],1));
r = unique(reshape(table2array(motifs(:,rxnRoles)),[],1));
metNames = string(model.mets(:)); rxnNames = string(model.rxns(:));
if isfield(model,'metNames'), metNames = string(model.metNames(:)); end
if isfield(model,'rxnNames'), rxnNames = string(model.rxnNames(:)); end
formulas = strings(numel(model.mets),1);
if isfield(model,'metFormulas') && numel(model.metFormulas) == numel(formulas)
    formulas = string(model.metFormulas(:));
end
gpr = strings(numel(model.rxns),1);
if isfield(model,'grRules'), gpr = string(model.grRules(:)); end
equations = strings(numel(r),1);
for k = 1:numel(r)
    substrates = find(model.S(:,r(k)) < 0);
    products = find(model.S(:,r(k)) > 0);
    equations(k) = localSide(substrates,model.S(:,r(k)),string(model.mets(:))) + ...
        " -> " + localSide(products,model.S(:,r(k)),string(model.mets(:)));
end
checkpoint.Differentiator_Metabolites = table(m,string(model.mets(m)),metNames(m),formulas(m), ...
    'VariableNames',{'Index','ID','Name','Formula'});
checkpoint.Differentiator_Reactions = table(r,string(model.rxns(r)),rxnNames(r),codes(r),gpr(r),equations, ...
    'VariableNames',{'Index','ID','Name','EC','GPR','ForwardEquation'});
% Tables/compatibility structs and model.S are deliberately not duplicated.
end

function side = localSide(indices,coefficients,ids)
if isempty(indices), side = "empty"; return; end
terms = strings(numel(indices),1);
for k = 1:numel(indices)
    coefficient = abs(full(coefficients(indices(k))));
    if coefficient == 1
        terms(k) = ids(indices(k));
    else
        terms(k) = string(coefficient)+" "+ids(indices(k));
    end
end
side = strjoin(terms," + ");
end
