% Write a function to pull out essential genes from the single gene
% deletion screen

% Take index of rows in grRateKO with value of 0 and output the
% corresponding genes as essential

% Function: essentialGenes(modelName, grRateKO)
    % Inputs:
        % modelName: name of your FBA model
        % grRateKO: Deletion strain growth rates (1/h) 
        % grThreshold: growth rate (1/h) cutoff for a complete KO
    % Outputs:
        % EG = cell array of the names of essential genes
        % Essential genes have a grRateKO = 0
        
function [egNames, egNum] = essentialGenes(modelName, grRateKO, grThreshold)
% get the gene names from your model
allGenes = modelName.genes;

% Create an empty cell array, EG
eg = {}; 
% start count to index how genes are added to EG
EG_count = 0;

% identify the indexes of zeros in grRateKO
% add gene(index) to EG
for i = 1:length(grRateKO);
    if grRateKO(i) <= grThreshold;
        EG_count = EG_count + 1;
        eg{EG_count} = allGenes(i);
    end 
end 
egNames = eg';
egNum = EG_count;



