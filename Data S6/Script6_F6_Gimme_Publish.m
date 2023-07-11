% Script to prepare a GIMME model using gene expression (RNA Seq) data and
% proteomics (iBAQ) data
% Written by Isabella Casini
% Last Updated: 2022.07.20

% USER INPUT: == there is a different line per microbes to be changed by the user

% USER INPUT: Only need initCobraToolbox to run once
% initCobraToolbox

solverOK = changeCobraSolver('gurobi', 'all', 1);
solverVersion = getCobraSolverVersion('gurobi', 1, 'C:\gurobi901');

% cd into the model directory
cd F:\GEM_Models\COBRAToolbox\F6_Gimme_2

% reads the model
% USER INPUT: Change depending on the microbe
% model = readCbModel('20220623_DH_model_bound_adjval.mat');
% model = readCbModel('20220623_ZZ_model_bound_adjval.mat');
model = readCbModel('20220623_MM_model_bound_adjval.mat');

changeCobraSolver ('gurobi', 'all'); % glpk, ibm_cplex, gurobi

% set the objective function to biomass, if not already
% model = changeObjective(model,'EX_biomass_e',1);                                    %ATP hydrolysis reaction 

% constrain the NGAM reaction, if not already
% model = changeRxnBounds(model,'ATPM',1.5,'l');
% model = changeRxnBounds(model,'ATPM',1000,'u');

% Make FRH irreversible for growth on H2/CO2 (equilibrator), if not already
model = changeRxnBounds(model,'FRH',0,'l');
model = changeRxnBounds(model,'FRH',1000,'u');

% Make Eha/Ehb irreversible for growth on H2/CO2 (equilibrator), if not
% already
% model = changeRxnBounds(model,'Eha',0,'l');
% model = changeRxnBounds(model,'Eha',1000,'u');
% model = changeRxnBounds(model,'Ehb',0,'l');
% model = changeRxnBounds(model,'Ehb',1000,'u');

% %----------------------------------------------------------------
% % PREPROCESSING DATA
% % https://opencobra.github.io/cobratoolbox/stable/modules/dataIntegration/transcriptomics/preprocessing/index.html
% 
% % defining gene expression data
% % first row: each column is a gene
% % second row: each column has the respective gene's expression data
% 
% % define path to gene data 
% USER INPUT: CHANGE FOR EACH MICROBE
% genedata_name = "F:/Experimental_Work/Transcriptomics/F6_Exp/salmon_DH/F6_RNA_DH_AVG_matlab.txt";
% genedata_name = "F:/Experimental_Work/Transcriptomics/F6_Exp/Updated_salmon_ZZ/F6_RNA_ZZ_AVG_matlab.txt";
genedata_name = "F:/Experimental_Work/Transcriptomics/F6_Exp/Updated_salmon_MM/F6_RNA_MM_AVG_matlab.txt";


% % loads the transcriptomics data
genedata = importdata(genedata_name);
genedata2 = genedata; % makes copy to edit
% genedata will end up with 3 structures
% data: [1×523 double]
% textdata: {1×523 cell}
%colheaders: {1×523 cell}

% Change the structure to have the row with genes called "gene"
genedata2 = renameStructField(genedata2,'textdata','gene');

% Change the structure to have the row with the data called "value"
genedata2 = renameStructField(genedata2,'data','value');

[expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model, genedata2,'True');

% USER INPUT: CHANGE FOR EACH MICROBE
% protdata_name = "F:/Experimental_Work/Proteomics/F6_2_results/Proteome_Discoverer/F6_Prot_DH_IBAQ_AVG_matlab.txt";
% protdata_name = "F:/Experimental_Work/Proteomics/F6_2_results/Proteome_Discoverer/F6_Prot_ZZ_IBAQ_AVG_matlab.txt";
protdata_name = "F:/Experimental_Work/Proteomics/F6_2_results/Proteome_Discoverer/F6_Prot_MM_IBAQ_AVG_matlab.txt";

% loads the proteomics data
protdata = importdata(protdata_name);
protdata2 = protdata; % makes copy to edit
% genedata will end up with 3 structures
% data: [1×523 double]
% textdata: {1×523 cell}
%colheaders: {1×523 cell}

% Change the structure to have the row with genes called "gene"
protdata2 = renameStructField(protdata2,'textdata','gene');

% Change the structure to have the row with the data called "value"
protdata2 = renameStructField(protdata2,'data','value');

[expressionRxnsprot, parsedGPRprot, gene_usedprot] = mapExpressionToReactions(model, protdata2,'True');

%---------------------------------------------------------------------
% GIMME
% set threshold

% cutoff at lower quartile
% cutoff  = quantile(gene_exp(~isnan(gene_exp)),0.25);
threshold_T  = quantile(expressionRxns(~isnan(expressionRxns)),0.25)
threshold_P  = quantile(expressionRxnsprot(~isnan(expressionRxnsprot)),0.25)

obj_frac = 0.9;

model_gimme_trans = GIMME(model, expressionRxns, threshold_T, obj_frac)
model_gimme_prot = GIMME(model, expressionRxnsprot, threshold_P, obj_frac)

% % %---------------------------------------------------------------------

% to run FBA
solution = optimizeCbModel(model)
solution_trans = optimizeCbModel(model_gimme_trans)
solution_prot = optimizeCbModel(model_gimme_prot)

% to run FBA without loops
solution_LL = optimizeCbModel(model,'max',0,0)
solution_transLL = optimizeCbModel(model_gimme_trans,'max',0,0)
solution_protLL = optimizeCbModel(model_gimme_prot,'max',0,0)

% Try with the objective function to maximize ATPM
model =  changeObjective(model,'ATPM',1);
model_gimme_trans = changeObjective(model_gimme_trans,'ATPM',1);                                    %ATP hydrolysis reaction 
model_gimme_prot = changeObjective(model_gimme_prot,'ATPM',1);                                    %ATP hydrolysis reaction 

% to run FBA
solution_ATP = optimizeCbModel(model)
solution_transATP = optimizeCbModel(model_gimme_trans)
solution_protATP = optimizeCbModel(model_gimme_prot)

% to run FBA without loops
solution_LLATP = optimizeCbModel(model)
solution_transLLATP = optimizeCbModel(model_gimme_trans,'max',0,0)
solution_protLLATP = optimizeCbModel(model_gimme_prot,'max',0,0)

% write out the data in a structured way
% https://www.mathworks.com/matlabcentral/answers/341935-combine-doubles-and-categorical

% define the outpath
% USER INPUT: Change path (and output filename)
filename = sprintf('F:/Experimental_Work/Reactors_Ley/F6_Exp2/FBA6/GIMME_results_MM_%s.xlsx', datestr(now,'yyyy-mm-dd_HH-MM'))

% normal model
rownames = {'rxn';'bio';'ATP';'bio_LL';'ATP_LL'};
T = table(model.rxns,solution.v,solution_ATP.v,solution_LL.v,solution_LLATP.v,'VariableNames',rownames);

% Transcriptomics model
rownamesT = {'rxn';'Trans_bio';'Trans_ATP';'Trans_bio_LL';'Trans_ATP_LL'};
TT = table(model_gimme_trans.rxns,solution_trans.v,solution_transATP.v,solution_transLL.v,solution_transLLATP.v,'VariableNames',rownamesT);

% Proteomics model
rownamesP = {'rxn';'Prot_bio';'Prot_ATP';'Prot_bio_LL';'Prot_ATP_LL'};
TP = table(model_gimme_prot.rxns,solution_prot.v,solution_protATP.v,solution_protLL.v,solution_protLLATP.v,'VariableNames',rownamesP);

% USER INPUT: CHANGE FOR EACH MICROBE

% sheet = 'DH_model';
% writetable(T,filename,'sheet',sheet,'Range','A1')
% sheet = 'DH_model_trans';
% writetable(TT,filename,'sheet',sheet,'Range','A1')
% sheet = 'DH_model_prot';
% writetable(TP,filename,'sheet',sheet,'Range','A1')

% sheet = 'ZZ_model';
% writetable(T,filename,'sheet',sheet,'Range','A1')
% sheet = 'ZZ_model_trans';
% writetable(TT,filename,'sheet',sheet,'Range','A1')
% sheet = 'ZZ_model_prot';
% writetable(TP,filename,'sheet',sheet,'Range','A1')

sheet = 'MM_model';
writetable(T,filename,'sheet',sheet,'Range','A1')
sheet = 'MM_model_trans';
writetable(TT,filename,'sheet',sheet,'Range','A1')
sheet = 'MM_model_prot';
writetable(TP,filename,'sheet',sheet,'Range','A1')
