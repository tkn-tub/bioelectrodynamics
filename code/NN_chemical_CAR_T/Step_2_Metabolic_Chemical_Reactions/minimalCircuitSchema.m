function s = minimalCircuitSchema(stage,mapCsLeak)
% Current 5+5 differentiators. Regulatory arrows do not consume metabolites.
% Cs support is external; mapCsLeak=true additionally maps basal Cs removal.
if nargin < 2, mapCsLeak = false; end
validateattributes(stage,{'numeric'},{'scalar','integer','>=',2,'<=',4});
if ~isscalar(mapCsLeak) || ~ismember(mapCsLeak,[0 1])
    error('MinimalCircuit:BufferMode','MAP_COUNTER_CS_LEAK must be true or false.');
end
s.mets = {'m_in','m_M1','m_Y1','m_F1_delay','m_M2','m_Y2'};
s.rxns = {'r_S1','r_S2','r_F','r_I','r_R','r_S1_2','r_S2_2','r_F2','r_I2','r_R2'};
s.consumes = logical([1 0 1 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0 0 0; ...
    0 0 0 1 0 1 0 1 0 0;0 0 0 1 1 0 0 0 0 0; ...
    0 0 0 0 0 0 1 0 1 0;0 0 0 0 0 0 0 0 1 1]);
s.produces = logical([0 0 0 0 0 0 0 0 0 0;1 0 0 0 0 0 0 0 0 0; ...
    0 1 0 0 0 0 0 0 0 0;0 0 1 0 0 0 0 0 0 0; ...
    0 0 0 0 0 1 0 0 0 0;0 0 0 0 0 0 1 1 0 0]);
s.regulators = {}; s.targets = {}; s.effects = zeros(1,0);
if stage >= 3
    s.mets = [s.mets,{'m_c','m_Cs'}];
    counterRxns = {'r_Cp','r_Cg','r_Cr'};
    counterConsumes = logical([0 1 0;0 1 1]);
    counterProduces = logical([1 0 0;0 0 0]);
    if mapCsLeak
        counterRxns{end+1} = 'r_CsLeak';
        counterConsumes(:,end+1) = [false;true];
        counterProduces(:,end+1) = false;
    end
    s.rxns = [s.rxns,counterRxns];
    s.consumes = blkdiag(s.consumes,counterConsumes);
    s.produces = blkdiag(s.produces,counterProduces);
    s.regulators = {'m_d2','m_d2'}; % m_d2=Y2, sole counter-control signal
    s.targets = {'r_Cp','r_Cr'};
    s.effects = [-1 -1];
end
if stage >= 4
    s.mets = [s.mets,{'m_Ni','m_Na','m_Nb','m_out'}];
    s.rxns = [s.rxns,{'r_Ni','r_NiDecay','r_Na','r_NaDecay', ...
        'r_Nb','r_NbDecay','r_No','r_NoDecay'}];
    s.consumes = blkdiag(s.consumes,logical([0 1 0 0 0 0 0 0; ...
        0 0 0 1 0 0 0 0;0 0 0 0 0 1 0 0;0 0 0 0 0 0 0 1]));
    s.produces = blkdiag(s.produces,logical([1 0 0 0 0 0 0 0; ...
        0 0 1 0 0 0 0 0;0 0 0 0 1 0 0 0;0 0 0 0 0 0 1 0]));
    s.regulators = [s.regulators,{'m_c','m_Ni','m_Ni','m_Na','m_Nb'}];
    s.targets = [s.targets,{'r_Ni','r_Na','r_Nb','r_No','r_No'}];
    s.effects = [s.effects,ones(1,5)];
end
s.dirs = strrep(s.rxns,'r_','dir_');
end
