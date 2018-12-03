% S3_Code_Implementation.m
% Run this to generate the data from S4 Data
% LJD, 12/12/17
% Goal: To perform gene essentiality simulations on the PA14 model on
% minimal media with carbon sources from Biolog Phenotypic Microarray Plates PM1 and PM2a. 
% Output: CSV file (titled '')of growth rate ratios for every gene on each
% carbon source. This data is equivalent to S4 Data. 

%% Initialize the workspace
clear all;
close all;
% Initialize the Cobra toolbox
initCobraToolbox;
% Set to the appropriate solver (see changeCobraSolver for solver options)
changeCobraSolver('gurobi5');
% Load the model, iPau1129 (Bartell, Blazier, et al 2017)
load('model_PA.mat');

%% Perform gene essentiality screens for compounds with exchange reactions in the model:

% Compound IDs (model formality)
limEX = {'EX_cpd00029(e)'...
'EX_cpd00137(e)'...
'EX_cpd00080(e)'...
'EX_cpd00117(e)'...
'EX_cpd00082(e)'...
'EX_cpd00182(e)'...
'EX_cpd00106(e)'...
'EX_cpd00051(e)'...
'EX_cpd00132(e)'...
'EX_cpd00041(e)'...
'EX_cpd00023(e)'...
'EX_cpd00053(e)'...
'EX_cpd00119(e)'...
'EX_cpd00100(e)'...
'EX_cpd00033(e)'...
'EX_cpd00380(e)'...
'EX_cpd00064(e)'...
'EX_cpd00066(e)'...
'EX_cpd00129(e)'...
'EX_cpd00054(e)'...
'EX_cpd00308(e)'...
'EX_cpd00035(e)'...
'EX_cpd00027(e)'...
'EX_cpd00118(e)'...
'EX_cpd00020(e)'...
'EX_cpd00036(e)'...
'EX_cpd00024(e)'...
'EX_cpd00281(e)'...
'EX_cpd00322(e)'...
'EX_cpd00130(e)'...
'EX_cpd00159(e)'...
'EX_cpd00107(e)'...
'EX_cpd00314(e)'...
'EX_cpd00162(e)'...
'EX_cpd00072(e)'...
'EX_cpd00280(e)'...
'EX_cpd00089(e)'...
'EX_cpd00079(e)'...
'EX_cpd00609(e)'...
'EX_cpd00164(e)'...
'EX_cpd00588(e)'...
'EX_cpd00154(e)'...
'EX_cpd00139(e)'...
'EX_cpd11589(e)'...
'EX_cpd00060(e)'...
'EX_cpd01242(e)'...
'EX_cpd00184(e)'...
'EX_cpd00156(e)'...
'EX_cpd00179(e)'...
'EX_cpd01262(e)'...
'EX_cpd00121(e)'...
'EX_cpd00652(e)'...
'EX_cpd00224(e)'...
'EX_cpd00142(e)'...
'EX_cpd00386(e)'...
'EX_cpd00105(e)'...
'EX_cpd00550(e)'...
'EX_cpd00246(e)'...
'EX_cpd00161(e)'...
'EX_cpd11592(e)'...
'EX_cpd11588(e)'...
'EX_cpd00851(e)'...
'EX_cpd00211(e)'...
'EX_cpd11585(e)'...
'EX_cpd00039(e)'...
'EX_cpd00489(e)'...
'EX_cpd00141(e)'...
'EX_cpd00249(e)'...
'EX_cpd00797(e)'...
'rJB00280'...
'EX_cpd01307(e)'...
'EX_cpd00122(e)'...
'EX_cpd00108(e)'...
'EX_cJB00034(e)'...
'EX_cpd00306(e)'}; 

% Compound names (Consistent with Phenotypic Microarrays)
limEX_Names = {'Acetic Acid'...
'Citric Acid'...
'D,L-alpha-Glycerol Phosphate'...
'D-Alanine'...
'D-Fructose'...
'Adenosine'...
'Fumaric Acid'...
'L-Arginine'...
'L-Asparagine'...
'L-Aspartic Acid'...
'L-Glutamic Acid'...
'L-Glutamine'...
'L-Histidine'...
'Glycerol'...
'Glycine'...
'Itaconic Acid'...
'L-Ornithine'...
'L-Phenylalanine'...
'L-Proline'...
'L-Serine'...
'Malonic Acid'...
'L-Alanine'...
'alpha-D-Glucose'...
'Putrescine'...
'Pyruvic Acid'...
'Succinic Acid'...
'alpha-Keto-Glutaric Acid'...
'gamma-Amino Butyric Acid'...
'L-Isoleucine'...
'L-Malic Acid'...
'L-Lactic Acid'...
'L-Leucine'...
'D-Mannitol'...
'2-Aminoethanol'...
'D-Fructose-6-Phosphate'...
'D-Galacturonic Acid'...
'D-Glucose-1-Phosphate'...
'D-Glucose-6-Phosphate'...
'D-Saccharic Acid'...
'D-Glucuronic Acid'...
'D-Sorbitol'...
'D-Xylose'...
'Glycolic Acid'...
'Glycyl-L-Aspartic Acid'...
'L-Methionine'...
'2-Deoxy-D-Ribose'...
'Thymidine'...
'L-Valine'...
'Maltose'...
'Maltotriose'...
'm-Inositol'...
'Mucic Acid'...
'L-Arabinose'...
'Acetoacetic Acid'...
'D-Malic Acid'...
'D-Ribose'...
'D-Serine'...
'Inosine'...
'L-Threonine'...
'Glycyl-L-Glutamic Acid'...
'Glycyl-L-Proline'...
'Hydroxy-L-Proline'...
'Butyric Acid'...
'L-Alanyl-Glycine'...
'L-Lysine'...
'p-Hydroxy Phenyl Acetic Acid'...
'Propionic Acid'...
'Uridine'...
'beta-Hydroxy Butyric Acid'...
'D-Gluconic Acid'...
'D-Arabitol'...
'N-Acetyl-D Glucosamine'...
'D-Galactose'...
'L-Sorbose'...
'Xylitol'}; 

% Initialize matrices to collect results
sol_MM = zeros(1,length(limEX));
grRatios = [];

% Collect results
for i = 1:length(limEX)
    % Set the media condition to MM + carbon source
    model_MM = changeMinimalMedia(model_PA, {limEX{i}});
    % Test that the base model can grow on each carbon source:
    sol = optimizeCbModel(model_MM);
    sol_MM(i) = sol.f;
    % Perform gene essentiality and save the growth ratios (grRatio)
    grRatio = singleGeneDeletion(model_MM);
    grRatios = [grRatios, grRatio];
end 

%% Add exchange reactions for compounds that are in biolog and model but don't have exchanges but have extracellular metabolite
% I.e. these compounds have transport reactions but no exchange reactions.
% Exchange reactions were temporarily added 
% Assumption: Compounds lacking a transporter were assumed not to support
% growth.

% Compound IDs
metList = {'cpd00136[e]'...
    'cpd00266[e]'...
    'cpd00138[e]'...
    'cpd00666[e]'...
    'cpd00047[e]'...
    'cpd00666[e]'...
    'cpd00794[e]'...
    'cpd00666[e]'...
    'cpd01502[e]'};
lb = ones(1,length(metList)).*-1000; 
ub = ones(1,length(metList)).*1000;

% Compound names
metNames = {'4-Hydroxy Benzoic Acid'...
    'D,L-Carnitine'...
    'D-Mannose'...
    'D-Tartaric Acid'...
    'Formic Acid'...
    'L-Tartaric Acid'...
    'D-Trehalose'...
    'm-Tartaric Acid'...
    'Citraconic Acid'};

% For all of the metabolites:
    % Add an exchange reaction
    % Change media to minimal media + that exchange reaction
    % Perform gene essentiality

% Initialize matrix to add growth ratios to
sol_MM2 = zeros(1,length(metList));
grRatios_EX = [];

% Collect results
for j = 1:length(metList)
    % Add exchange reaction
    [model_addEX, name_EX] = addExchangeRxn_JB(model_PA, {metList{j}}, [lb(j)], [ub(j)]);
    % Change media to the new exchange reaction
    model_MM = changeMinimalMedia(model_addEX,name_EX);
    % Test that the base model can grow on each carbon source:
    sol = optimizeCbModel(model_MM);
    sol_MM2(j) = sol.f;
    % Perform gene essentiality and save the essential genes
    grRatio = singleGeneDeletion(model_MM);
    grRatios_EX = [grRatios_EX, grRatio];
end 

%% Export all predictions 
varNames = {'V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V13','V14','V15','V16','V17','V18','V19','V20','V21','V22','V23','V24','V25','V26','V27','V28','V29','V30','V31','V32','V33','V34','V35','V36','V37','V38','V39','V40','V41','V42','V43','V44','V45','V46','V47','V48','V49','V50','V51','V52','V53','V54','V55','V56','V57','V58','V59','V60','V61','V62','V63','V64','V65','V66','V67','V68','V69','V70','V71','V72','V73','V74','V75','V76','V77','V78','V79','V80','V81','V82','V83','V84','V85'};
headings = [limEX_Names(1,:), metNames(1,:)];
grRatios_all = [grRatios, grRatios_EX];

% Make a table with all information 
tableData = table(model_PA.genes, 'VariableNames',{'V0'});
tableGenes = table({'gene'}, 'VariableNames',{'V0'});
for i = 1:84
    a = table(model_PA.genes,num2cell(grRatios_all(:,i)), 'VariableNames',{'V0',varNames{i}});
    b = table({'gene'}, headings(i), 'VariableNames', {'V0',varNames{i}});
    tableData = join(tableData, a);
    tableGenes = join(tableGenes, b);
end 
tableComplete = vertcat(tableGenes, tableData);

% Export to a CSV file
writetable(tableComplete, 'geneEssentialityPredictions.csv','Delimiter', ',', 'QuoteStrings', true, 'WriteVariableNames',false);
