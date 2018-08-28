% Window Size Analysis
% LJD
% 8/23/18
% Goal: Evaluate the effect of changing the window size on growth rate
% calculations 
close all;
clear all;
load('S1_RawData.mat');

%% Looking at window size and doubling time
% Generate sample data
x = [1:1:272].*10./60;
n = [4,2, 1.15, 1.07]'; % corresponding to DT of 30min, 1hr,5hr, 10hr
y_true = n.^x;
numCurves = length(n);
windows = 2:2:64;
numTrials = 100;

runningGrowthRate_true = zeros(numCurves,length(windows));
runningGrowthRate_measured = zeros(numCurves,length(windows));
for numRandNoise = 1:1:numTrials
y = n.^x + ((0.1.*n.^x).*randi([-1,1],1,length(x)));
growthRates2 = []; 
growthRates3 = [];
% Iterate over many window sizes. What is the DT?
    for i = windows
        sampleData2 = struct('time',x,'data',y,'wellNames',{wellNames(1:numCurves)});
        sampleData3 = struct('time',x,'data',y_true,'wellNames',{wellNames(1:numCurves)});
        growthDynamStruct2 = calcGrowthDynamics(sampleData2, i);
        growthDynamStruct3 = calcGrowthDynamics(sampleData3, i);
        growthRates2 = [growthRates2, growthDynamStruct2.growthRate];
        growthRates3 = [growthRates3, growthDynamStruct3.growthRate];
    end 
 runningGrowthRate_true = runningGrowthRate_true + growthRates3;  
 runningGrowthRate_measured = runningGrowthRate_measured + growthRates2;  

end 
avgGR_true = runningGrowthRate_true./numTrials;
avgGR_measured = runningGrowthRate_measured./numTrials;
measured_dt_true = log(2)./avgGR_true;
measured_dt = log(2)./avgGR_measured;

percent_error = abs(measured_dt - measured_dt_true)./measured_dt_true.*100;
%% Plot the results
close all;
% For window sizes 30min, 5h, 10h, 20h, 35h
% load('windowSize_100RandomTrials.mat');
% For window sizes 30min, 1h, 5h, 10h
load('windowSize_100RandomTrials_10h.mat');
titles = {'DT = 30m','DT = 1h','DT = 5h','DT = 10h'};
figure();hold on;

% Plot the DT measured and the true doubling time by window size
% Also plot the percent error between the two for each window size
ncol = length(n);
for i = 1:length(n)
%     subplot(2,ncol,i); hold on;
%     scatter([2:2:64].*10./60, measured_dt(i,:));hold on; plot([2:2:64].*10./60, measured_dt_true(i,:));xlabel('window size (hr)');
% ylabel('doubling time (hr)');%legend('doubling time 10% noise','true doubling time');
% 
    
    subplot(1,ncol,i); hold on;
    scatter([2:2:64].*10, percent_error(i,:)); xlabel('window size (min)');
    ylabel('percent error abs(calc - actual)/actual * 100');ylim([0,50]);hold on;
   title(titles{i})
end 

% Plot error for selected window size 280 minutes (28 time points)
figure();
plot(measured_dt_true(1:ncol,14), percent_error(:,14));hold on;xlabel('Doubling Time');
    ylabel('Percent Error abs(calc - actual)/actual * 100');ylim([0,10]);
    legend(cellstr(num2str((windows(14)).*10, 'DT=%-d')));title('Error vs Doubling time for Window Size 280 minutes');
% 
%     % Plot DT vs percent error for different window sizes
% figure();
% for j = 1:length(windows)
%     plot(measured_dt_true(1:ncol,j), percent_error(:,j));hold on;xlabel('Doubling Time');
%     ylabel('percent error');
% end 
% hold on; scatter(measured_dt_true(1:ncol,14), percent_error(:,14),'filled')
% legendCell = cellstr(num2str((windows)', 'N=%-d'));
% legend(legendCell)
