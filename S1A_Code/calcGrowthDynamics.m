function growthDynamicsStruct = calcGrowthDynamics(growthStruct,generateFigures)
% S1 Code
% Dunphy, Yen, and Papin 2018
% Calculate the growth rate and time to mid-exponential phase from raw
% growth curve data

% Method adapted from: 	
% Hall BG, Acar H, Nandipati A, Barlow M. Growth Rates Made Easy. 
    % Mol Biol Evol. 2014 Jan 1;31(1):232–8. 

% Inputs: 
    % growthStruct - structure with the following fields:
        % data - n x m matrix where 
            % n = number of wells measured
            % m = number of timepoints
            % values in the matrix are optical density (OD) measurements
        % time - 1 x m matrix with timepoints for each measurement
        % wellNames - n x 1 matrix with labels for each well. 
            % The example code well names are carbon sources from Biolog
            % Phenotypic Microarray plates
    % generateFigures - binary call to plot individual curves and growth
    % rate locations (1) or omit plots (0). Default = 0. 

% Output: GRbiologStruct ->  structure with growth rates for each well
% added under a field, "growthRate", and time to mid-exponential added under field
% "lagPhase"

if nargin > 1
  genFigures = generateFigures;
else
  genFigures = 0;
end

% Set up growth rate and lag phase parameters to be filled in
growthStruct.growthRate = zeros(length(growthStruct.wellNames),1);
growthStruct.lagPhase = zeros(length(growthStruct.wellNames),1);

% Calculate the natural log (ln) of the raw OD600 data
growthStruct.lnData = log(growthStruct.data);

% Sliding window: calculate slope in windowSize point intervals
for j = 1:length(growthStruct.wellNames)    
    % Initialize search space
    slopes = []; % will include all slopes calculated across the growth curves
    slopesFiltered = []; % will only include slopes calculated prior to stationary phase
    linRange = []; % will include the indices of windows in exponential phase
    
    % Set parameters
    % These parameters were optimized by comparing to growth rate measures
    % by hand
    windowSize = 8; % number of points included in each sliding window calculation
    slopeThresh = 0.75; % Multiplied by the max growth rate to determine other exponential phase windows
    neighborThresh = 15; % max distance apart in linrange points can be
    stationaryThresh = 1; % fraction of max lnOD after which we assume stationary
    
    % Identify valid time range to identify exponential phase by finding
    % ~start of stationary phase. We will calculate the growth rate before
    % stationary phase
    stationaryTime = find(growthStruct.lnData(j,:) >= log(stationaryThresh*max(growthStruct.data(j,:))));
    if stationaryTime(1) > 1 % if max OD is not at the first timepoint
        stationaryTime = stationaryTime(1);
    else % recalulate stationary time until get past at least one point not above max...not sure about this part...
        index = 1;
        while stationaryTime(1) == 1;
            if index < length(growthStruct.time)
                index = index + 1;
                stationaryTime = find(growthStruct.lnData(j,index:end) >= log(stationaryThresh*max(growthStruct.data(j,index:end))));
            else % if you get to the last point without finding a max (happens if curve is very flat or even decreasing), just don't identify stationary and use whole space. 
                stationaryTime = [length(growthStruct.time)];
            end
        end
        stationaryTime = stationaryTime(1);
    end 

    % Implement the sliding window algorithm
    % Calculates slope of time vs ln(OD600) in windowSize increments across
    % the growth curve
    for i = 1:length(growthStruct.time)-(windowSize-1)
        time = growthStruct.time(i:i+(windowSize-1));
        data = growthStruct.lnData(j,i:i+(windowSize-1));
        m = polyfit(time, data, 1);
        % Save slope
        slopes = [slopes, m(1)];
    end 
    
    % Filter slopes so that we only consider those prior to stationary
    % phase
    if stationaryTime < length(slopes) 
        slopesFiltered = slopes(1:stationaryTime);
    else
        slopesFiltered = slopes; % prevents code from breaking if last point denotes start of stationary
    end
    
    % Identify windows where slope was >= slopeThresh of max slope to find
    % the exponential range
    % If slopes are all negative, linRange will be empty.
    linRange = slopesFiltered >= slopeThresh*max(slopesFiltered);
    linRange = find(double(linRange)); % indicies (start of windows) where slope >= 95% max slope
 
    % Limit linRange to neighboring (within neighborThresh) points:
    closeness = diff(linRange);
    if max(closeness)>neighborThresh
        locNeighbors = find(closeness>neighborThresh);
        rangeLengths = zeros(1,length(locNeighbors)+1);
        ranges = zeros(length(locNeighbors)+1,2);
        for q = 1:length(locNeighbors)
            if q == 1
                rangeLength = length(linRange(1:locNeighbors(q)));
                ranges(q,:) = [1,locNeighbors(1)];
            else
                rangeLength = length(linRange(locNeighbors(q-1)+1:locNeighbors(q)));
                ranges(q,:) = [locNeighbors(q-1)+1, locNeighbors(q)];
            end 
            rangeLengths(q) = rangeLength;
        end 
        rangeEnd = length(linRange(locNeighbors(end)+1:end));
        rangeLengths(end) = rangeEnd;
        ranges(end,:) = [locNeighbors(end)+1,length(linRange)];
        % Use first instance of a range of maximum length
        a = find(rangeLengths == max(rangeLengths));
        ranges(a(1),:);
        linRange = linRange(ranges(a(1),:));
    end

    % Calculate the growth rate of the exponential range (using slope of
    % polyfit) and the start of mid exponential phase (first timepoint of
    % exponential range)
    if isempty(linRange) % catches issue of linRange being empty due to all slopes before stationary being < 0
        warning('Problem finding exponential phase. Assigning a lag phase time of 0 and a growth rate of 0')
        errorLocation = {growthStruct.wellNames{j}}
        % Assume in this case that there was no growth
        growthStruct.growthRate(j,1) = 0;
        growthStruct.lagPhase(j,1) = 0;
    else
        linRangeStart = min(linRange);
        linRangeStop = max(linRange) + 4;
        P = polyfit(growthStruct.time(linRangeStart:linRangeStop), growthStruct.lnData(j,linRangeStart:linRangeStop),1);
        growthStruct.growthRate(j,1) = P(1);
        growthStruct.lagPhase(j,1) = growthStruct.time(linRangeStart);
    end 
  
    %% uncomment if you want to plot each well...slow so i commented it out
   if genFigures == 1
    figure()
    subplot(3,1,1)
    scatter(growthStruct.time, growthStruct.data(j,:),5,'k'); hold on; ylabel('OD'); axis([0 48 0 1]);
    legend(num2str(max(growthStruct.data(j,:))));title(growthStruct.wellNames(j));
    
    subplot(3,1,2);
    scatter(growthStruct.time, growthStruct.lnData(j,:),5,'k'); hold on; ylabel('lnOD'); 
    plot(growthStruct.time,ones(1,length(growthStruct.time))*max(growthStruct.lnData(j,:))*stationaryThresh,'c'); scatter(growthStruct.time(stationaryTime), stationaryThresh*max(growthStruct.lnData(j,:)),'m');
    plot(growthStruct.time(linRangeStart:linRangeStop), growthStruct.time(linRangeStart:linRangeStop)*P(1) + P(2),'-r','LineWidth',5)
    axis([0 48 -2.5 0.1]); legend(num2str(max(growthStruct.lnData(j,:))));
    
    subplot(3,1,3);
    scatter(growthStruct.time(1:end-(windowSize-1)),slopes,5,'k'); 
    hold on;ylabel('dLnOD/dt');
   end 
   
end 
growthDynamicsStruct = growthStruct;

end 

