% S3_Code_For_Reviewers
% LJD 8/14/18
% Test the impact of varying the lb on the exchange reactions of 4-HBA,
% L-Leucine, and L-Isoleucine to see if the flux impacts essentiality
% predictions.

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
limEX = {'EX_cpd00107(e)','EX_cpd00322(e)'}; 

% Compound names (Consistent with Phenotypic Microarrays)
limEX_Names = {'L-Leucine','L-Isoleucine'}; 

% Initialize matrices to collect results
sol_MM_10 = zeros(1,length(limEX));
sol_MM_1000 = zeros(1,length(limEX));
grRatios_1000 = [];
grRatios_10 = [];

% Collect results
for i = 1:length(limEX)
    % Set the media condition to MM + carbon source
    % Make the lower bound -1000 instead of -10 to address reviewer
    % feedback
    model_MM_10 = changeMinimalMedia(model_PA, {limEX{i}},-10);
    model_MM_1000 = changeMinimalMedia(model_PA, {limEX{i}},-1000);

    % Test that the base model can grow on each carbon source:
    sol_10 = optimizeCbModel(model_MM_10);
    sol_1000 = optimizeCbModel(model_MM_1000);

    sol_MM_10(i) = sol_10.f;
    sol_MM_1000(i) = sol_1000.f;

    % Perform gene essentiality and save the growth ratios (grRatio)
    grRatio_10 = singleGeneDeletion(model_MM_10);
    grRatio_1000 = singleGeneDeletion(model_MM_1000);

    grRatios_10 = [grRatios_10, grRatio_10];
    grRatios_1000 = [grRatios_1000, grRatio_1000];
    [egNames_10, egNum_10] = essentialGenes(model_MM_10, grRatio_10, 0.0001);
    [egNames_1000, egNum_1000] = essentialGenes(model_MM_1000, grRatio_1000, 0.0001);
    [egNum_10, egNum_1000]
    isequal(egNames_1000, egNames_10)

end 

%% Add exchange reactions for compounds that are in biolog and model but don't have exchanges but have extracellular metabolite
% I.e. these compounds have transport reactions but no exchange reactions.
% Exchange reactions were temporarily added 
% Assumption: Compounds lacking a transporter were assumed not to support
% growth.

% Compound IDs
metList = {'cpd00136[e]'};
lb = ones(1,length(metList)).*-1000; 
ub = ones(1,length(metList)).*1000;

% Compound names
metNames = {'4-Hydroxy Benzoic Acid'};

% For all of the metabolites:
    % Add an exchange reaction
    % Change media to minimal media + that exchange reaction
    % Perform gene essentiality

% Initialize matrix to add growth ratios to
sol_MM2_10 = zeros(1,length(metList));
sol_MM2_1000 = zeros(1,length(metList));

grRatios_EX_10 = [];
grRatios_EX_1000 = [];


% Collect results
for j = 1:length(metList)
    % Add exchange reaction
    [model_addEX, name_EX] = addExchangeRxn_JB(model_PA, {metList{j}}, [lb(j)], [ub(j)]);
    % Change media to the new exchange reaction
    % Make the lower bound -1000 instead of -10 to address reviewer
    % feedback
    model_MM_10 = changeMinimalMedia(model_addEX,name_EX, -10);
    model_MM_1000 = changeMinimalMedia(model_addEX,name_EX, -1000);

    % Test that the base model can grow on each carbon source:
    sol_10 = optimizeCbModel(model_MM_10);
    sol_1000 = optimizeCbModel(model_MM_1000);

    sol_MM2_10(j) = sol_10.f;
    sol_MM2_1000(j) = sol_1000.f;

    % Perform gene essentiality and save the essential genes
    grRatio_10 = singleGeneDeletion(model_MM_10);
    grRatio_1000 = singleGeneDeletion(model_MM_1000);

    grRatios_EX_10 = [grRatios_EX_10, grRatio_10];
    grRatios_EX_1000 = [grRatios_EX_1000, grRatio_1000];
    [egNames_10, egNum_EX_10] = essentialGenes(model_MM_10, grRatio_10, 0.0001);
    [egNames_1000, egNum_EX_1000] = essentialGenes(model_MM_1000, grRatio_1000, 0.0001);
    [egNum_EX_10, egNum_EX_1000]
    isequal(egNames_1000, egNames_10)

end 
%% Compare an exchange of -10 (limited) vs and exchange of -1000 (open)
figure();
[n, xout] = hist([abs(grRatios_10 - grRatios_1000),abs(grRatios_EX_10 - grRatios_EX_1000)]);
bar(xout, n, 'barwidth', 1, 'basevalue', 0.1);
set(gca,'YScale','log');
xlabel('Relative change in gene effect (limited relative growth - open relative growth)'); ylabel('n');
legend('L-Leucine','L-Isoleucine','4-HBA');title('Bound of -10 to -1000');
