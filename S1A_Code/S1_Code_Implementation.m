% S1 Code Implementation - main script
% Dunphy, Yen, and Papin 2018
% Calculate and plot growth dynamics from raw growth curve data
% Lineages: Ancestor, LBE, PIP, CIP, TOB
% Carbon Sources: All 14 included in Figure 3

close all;
clear all;

% Load a sample dataset
load('S1_RawData.mat');
% This dataset structure (sampleData) contains the raw growth data all biological replicates
% of each lab-evolved lineage and the ancestral lineage. 
% Figure 5 and S3 Data contain values averaged across these three biological replicates. 

% Set up structure from raw data
sampleData = struct('time',time,'data',data,'wellNames',{wellNames});

% Calculate the growth dynamics (growth rate and time to mid exponential,
% 'lagPhase' for each growth curve)



%% Without figures
growthDynamics = calcGrowthDynamics(sampleData, 8);



% sampleGrowthDynamics is a new structure with all the info of sampleData
% plus growth rates and times to mid-exp for each growth curve.
% This structure was exported by the authors and the median was taken in R to generate
% S3 Data. 

% Uncomment below to see growth rates and times to mid exponential  
% growthDynamics.growthRate
% growthDynamics.lagPhase

% To output figures, uncomment the below line.
% Note: This may crash depending on graphics ability\
% Red line shows where slope was calculated. Left side of red line was
% determined to be the time to mid-exponential phase. 

% growthDynamics = calcGrowthDynamics(sampleData,8,1);


