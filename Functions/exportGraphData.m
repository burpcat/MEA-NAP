function exportGraphData(dataStage, expData, Params, ExN)
% Centralized function to export graph data at different pipeline stages
% CORRECTED VERSION - Fixed field names and added all required metrics
%
% Inputs:
%   dataStage - string indicating the pipeline stage ('spike', 'ephys', 'adjmatrix', 'netmet', 'nodecartography', 'all')
%   expData - structure containing experiment data (can include Info, adjMs, NetMet, etc.)
%   Params - pipeline parameters structure
%   ExN - experiment number (optional, only needed for certain stages)
%
% This function exports all data needed to recreate the 4B_GroupComparisons visualizations

% Check if graph data export is enabled
if ~isfield(Params, 'saveGraphData') || ~Params.saveGraphData
    return;
end

% Create graph data folder if it doesn't exist
if ~isfield(Params, 'graphDataFolder')
    Params.graphDataFolder = fullfile(Params.outputDataFolder, Params.outputDataFolderName, 'GraphData');
end

if ~isfolder(Params.graphDataFolder)
    mkdir(Params.graphDataFolder);
end

% Make sure we have the Info structure
if isfield(expData, 'Info')
    Info = expData.Info;
elseif isfield(expData, 'info')
    Info = expData.info;
else
    error('Info structure not found in expData');
end

% Get experiment name and group
if iscell(Info.FN)
    expName = char(Info.FN);
else
    expName = Info.FN;
end

if iscell(Info.Grp)
    groupName = char(Info.Grp);
else
    groupName = Info.Grp;
end

% Create folder structure
groupFolder = fullfile(Params.graphDataFolder, groupName);
if ~isfolder(groupFolder)
    mkdir(groupFolder);
end

expFolder = fullfile(groupFolder, expName);
if ~isfolder(expFolder)
    mkdir(expFolder);
end

% Export data based on stage
switch lower(dataStage)
    case 'spike'
        exportSpikeData(expData, expName, expFolder);
        
    case 'ephys'
        exportEphysData(expData, expName, expFolder, Params, ExN);
        
    case 'adjmatrix'
        exportAdjMatrixData(expData, expName, expFolder, Params);
        
    case 'netmet'
        exportNetMetData(expData, expName, expFolder, Params);
        
    case 'nodecartography'
        exportNodeCartographyData(expData, expName, expFolder, Params);
        
    case 'all'
        % Export everything available
        exportSpikeData(expData, expName, expFolder);
        exportEphysData(expData, expName, expFolder, Params, ExN);
        exportAdjMatrixData(expData, expName, expFolder, Params);
        exportNetMetData(expData, expName, expFolder, Params);
        exportNodeCartographyData(expData, expName, expFolder, Params);
        
    otherwise
        warning('Unknown data stage: %s', dataStage);
end

end

function exportSpikeData(expData, expName, saveFolder)
% Export spike detection data for raster plots and spike analysis

if isfield(expData, 'spikeTimes')
    spikeData = struct();
    spikeData.spikeTimes = expData.spikeTimes;
    
    if isfield(expData, 'channels')
        spikeData.channels = expData.channels;
    elseif isfield(expData, 'Info') && isfield(expData.Info, 'channels')
        spikeData.channels = expData.Info.channels;
    end
    
    if isfield(expData, 'Info')
        spikeData.Info = expData.Info;
    end
    
    filePath = fullfile(saveFolder, [expName, '_spikeData.mat']);
    save(filePath, 'spikeData');
    fprintf('Exported spike data: %s\n', filePath);
end

if isfield(expData, 'spikeMatrix')
    rasterData = struct();
    rasterData.spikeMatrix = expData.spikeMatrix;
    
    if isfield(expData, 'channels')
        rasterData.channels = expData.channels;
    elseif isfield(expData, 'Info') && isfield(expData.Info, 'channels')
        rasterData.channels = expData.Info.channels;
    end
    
    if isfield(expData.Info, 'duration_s')
        rasterData.duration = expData.Info.duration_s;
    end
    
    filePath = fullfile(saveFolder, [expName, '_spikeRaster.mat']);
    save(filePath, 'rasterData');
    fprintf('Exported spike raster data: %s\n', filePath);
end
end

function exportEphysData(expData, expName, saveFolder, Params, ExN)
% Export neuronal activity data - both electrode and recording level metrics
% Based on what PlotEphysStats.m actually uses

if ~isfield(expData, 'Ephys')
    return;
end

% ELECTRODE-LEVEL METRICS (for 2A_IndividualNeuronalAnalysis plots)
electrodeLevelData = struct();

% Corrected field names based on ExtractNetMet.m and pipeline code
electrodeFields = {
    'FR',                           % Firing rates
    'channelBurstRate',            % Channel burst rates  
    'channelBurstDur',             % Channel burst durations
    'channelFracSpikesInBursts',   % Fraction of spikes in bursts
    'channelISIwithinBurst',       % ISI within bursts
    'channeISIoutsideBurst'        % ISI outside bursts (note: typo in original code)
};

for i = 1:length(electrodeFields)
    fieldName = electrodeFields{i};
    if isfield(expData.Ephys, fieldName)
        electrodeLevelData.(fieldName) = expData.Ephys.(fieldName);
    else
        fprintf('Warning: Electrode field %s not found\n', fieldName);
    end
end

% Add spatial information
if isfield(expData, 'coords')
    electrodeLevelData.coords = expData.coords;
elseif isfield(Params, 'coords') && exist('ExN', 'var') && ExN <= length(Params.coords)
    electrodeLevelData.coords = Params.coords{ExN};
end

if isfield(expData, 'channels')
    electrodeLevelData.channels = expData.channels;
elseif isfield(expData.Info, 'channels')
    electrodeLevelData.channels = expData.Info.channels;
end

% Add experiment info
electrodeLevelData.Info = expData.Info;

% Save electrode-level data
filePath = fullfile(saveFolder, [expName, '_electrodeLevelActivity.mat']);
save(filePath, 'electrodeLevelData');
fprintf('Exported electrode-level activity data: %s\n', filePath);

% RECORDING-LEVEL METRICS (for PlotEphysStats group comparisons)
recordingLevelData = struct();

% All recording-level fields from PlotEphysStats.m NetMetricsE
recordingFields = {
    'numActiveElec',                    % Number of active electrodes
    'FRmean',                          % Mean firing rate
    'FRmedian',                        % Median firing rate  
    'FRstd',                           % Standard deviation of firing rate
    'FRsem',                           % Standard error of mean firing rate
    'FRiqr',                           % Interquartile range of firing rate
    'NBurstRate',                      % Network burst rate
    'meanNumChansInvolvedInNbursts',   % Mean channels in network bursts
    'meanNBstLengthS',                 % Mean network burst length (seconds)
    'meanISIWithinNbursts_ms',         % Mean ISI within network bursts
    'meanISIoutsideNbursts_ms',        % Mean ISI outside network bursts
    'CVofINBI',                        % CV of inter-network-burst intervals
    'fracInNburst',                    % Fraction of spikes in network bursts
    'channelAveBurstRate',             % Average channel burst rate
    'channelAveBurstDur',              % Average channel burst duration
    'channelAveISIwithinBurst',        % Average ISI within bursts
    'channelAveISIoutsideBurst',       % Average ISI outside bursts
    'channelAveFracSpikesInBursts',    % Average fraction spikes in bursts
    'numNbursts'                       % Number of network bursts
};

for i = 1:length(recordingFields)
    fieldName = recordingFields{i};
    if isfield(expData.Ephys, fieldName)
        recordingLevelData.(fieldName) = expData.Ephys.(fieldName);
    else
        fprintf('Warning: Recording field %s not found\n', fieldName);
    end
end

% Add experiment metadata
recordingLevelData.expName = expName;
recordingLevelData.Grp = expData.Info.Grp;
recordingLevelData.DIV = expData.Info.DIV;
if isfield(expData.Info, 'duration_s')
    recordingLevelData.duration_s = expData.Info.duration_s;
end

% Save recording-level data
filePath = fullfile(saveFolder, [expName, '_recordingLevelActivity.mat']);
save(filePath, 'recordingLevelData');
fprintf('Exported recording-level activity data: %s\n', filePath);

end

function exportAdjMatrixData(expData, expName, saveFolder, Params)
% Export adjacency matrices for all lag values

if ~isfield(expData, 'adjMs')
    return;
end

% Export for each lag value
for lagIdx = 1:length(Params.FuncConLagval)
    lagVal = Params.FuncConLagval(lagIdx);
    adjMatrixFieldName = sprintf('adjM%.fmslag', lagVal);
    
    if isfield(expData.adjMs, adjMatrixFieldName)
        adjMatrixData = struct();
        adjMatrixData.adjM = expData.adjMs.(adjMatrixFieldName);
        adjMatrixData.lagValue = lagVal;
        adjMatrixData.expName = expName;
        adjMatrixData.Info = expData.Info;
        
        % Add coordinates and channels if available
        if isfield(expData, 'coords')
            adjMatrixData.coords = expData.coords;
        end
        if isfield(expData, 'channels')
            adjMatrixData.channels = expData.channels;
        elseif isfield(expData.Info, 'channels')
            adjMatrixData.channels = expData.Info.channels;
        end
        
        filename = sprintf('%s_adjMatrix_lag%.f.mat', expName, lagVal);
        filePath = fullfile(saveFolder, filename);
        save(filePath, 'adjMatrixData');
        fprintf('Exported adjacency matrix (lag %.f): %s\n', lagVal, filePath);
    end
end
end

function exportNetMetData(expData, expName, saveFolder, Params)
% Export ALL network metrics - corrected field names based on ExtractNetMet.m
% This exports data needed for ALL 4B_GroupComparisons plots

if ~isfield(expData, 'NetMet')
    return;
end

% Export for each lag value
for lagIdx = 1:length(Params.FuncConLagval)
    lagVal = Params.FuncConLagval(lagIdx);
    lagFieldName = sprintf('adjM%.fmslag', lagVal);
    
    if ~isfield(expData.NetMet, lagFieldName)
        continue;
    end
    
    netMetLag = expData.NetMet.(lagFieldName);
    
    % NODE-LEVEL METRICS (for 1_NodeByGroup and 2_NodeByAge)
    nodeLevelData = struct();
    
    % Corrected field names from ExtractNetMet.m
    nodeFields = {
        'ND',           % Node degree
        'MEW',          % Mean edge weight  
        'NS',           % Node strength
        'Eloc',         % Local efficiency
        'Z',            % Within-module degree z-score
        'BC',           % Betweenness centrality
        'PC',           % Participation coefficient
        'aveControl',   % Average controllability
        'modalControl', % Modal controllability
        'NE',           % Nodal efficiency
        'activeChannel' % Active channel IDs
    };
    
    for i = 1:length(nodeFields)
        fieldName = nodeFields{i};
        if isfield(netMetLag, fieldName)
            nodeLevelData.(fieldName) = netMetLag.(fieldName);
        else
            fprintf('Warning: Node metric %s not found for lag %.f\n', fieldName, lagVal);
        end
    end
    
    % Add metadata
    nodeLevelData.lagValue = lagVal;
    nodeLevelData.expName = expName;
    nodeLevelData.Grp = expData.Info.Grp;
    nodeLevelData.DIV = expData.Info.DIV;
    
    % Add coordinates
    if isfield(expData, 'coords')
        nodeLevelData.coords = expData.coords;
    end
    if isfield(expData, 'channels')
        nodeLevelData.channels = expData.channels;
    elseif isfield(expData.Info, 'channels')
        nodeLevelData.channels = expData.Info.channels;
    end
    
    % Save node-level data
    filename = sprintf('%s_nodeLevelMetrics_lag%.f.mat', expName, lagVal);
    filePath = fullfile(saveFolder, filename);
    save(filePath, 'nodeLevelData');
    fprintf('Exported node-level metrics (lag %.f): %s\n', lagVal, filePath);
    
    % NETWORK-LEVEL METRICS (for 3_RecordingsByGroup and 4_RecordingsByAge)
    networkLevelData = struct();
    
    % All network-level metrics from ExtractNetMet.m
    networkFields = {
        'aN',               % Number of active nodes / network size
        'Dens',             % Density
        'NDmean',           % Mean node degree  
        'NDtop25',          % Top 25% node degree
        'sigEdgesMean',     % Significant edge weight mean
        'sigEdgesTop10',    % Top 10% edge weight mean
        'NSmean',           % Mean node strength
        'ElocMean',         % Mean local efficiency
        'CC',               % Clustering coefficient
        'nMod',             % Number of modules
        'Q',                % Modularity score
        'percentZscoreGreaterThanZero',  % Percentage within-module z-score > 0
        'percentZscoreLessThanZero',     % Percentage within-module z-score < 0  
        'PL',               % Mean path length
        'PCmean',           % Participation coefficient mean
        'PCmeanBottom10',   % Bottom 10% PC
        'PCmeanTop10',      % Top 10% PC
        'Eglob',            % Global efficiency
        'SW',               % Small-worldness sigma
        'SWw',              % Small-worldness omega  
        'aveControlMean',   % Mean average controllability
        'BCmeantop5',       % Top 5% betweenness centrality
        'Hub3',             % Hubs (3+ criteria)
        'Hub4'              % Hubs (4 criteria)
    };
    
    for i = 1:length(networkFields)
        fieldName = networkFields{i};
        if isfield(netMetLag, fieldName)
            networkLevelData.(fieldName) = netMetLag.(fieldName);
        else
            fprintf('Warning: Network metric %s not found for lag %.f\n', fieldName, lagVal);
        end
    end
    
    % Add lag-independent metrics (only for first lag)
    if lagIdx == 1
        lagIndependentFields = {
            'num_nnmf_components',      % Number NMF components
            'nComponentsRelNS',         % nNMF div network size  
            'effRank'                   % Effective rank
        };
        
        for i = 1:length(lagIndependentFields)
            fieldName = lagIndependentFields{i};
            if isfield(netMetLag, fieldName)
                networkLevelData.(fieldName) = netMetLag.(fieldName);
            end
        end
    end
    
    % Add metadata
    networkLevelData.lagValue = lagVal;
    networkLevelData.expName = expName;
    networkLevelData.Grp = expData.Info.Grp;
    networkLevelData.DIV = expData.Info.DIV;
    if isfield(expData.Info, 'duration_s')
        networkLevelData.duration_s = expData.Info.duration_s;
    end
    
    % Save network-level data
    filename = sprintf('%s_networkLevelMetrics_lag%.f.mat', expName, lagVal);
    filePath = fullfile(saveFolder, filename);
    save(filePath, 'networkLevelData');
    fprintf('Exported network-level metrics (lag %.f): %s\n', lagVal, filePath);
end
end

function exportNodeCartographyData(expData, expName, saveFolder, Params)
% Export node cartography data for 6_NodeCartographyByLag plots

if ~isfield(expData, 'NetMet')
    return;
end

% Export for each lag value
for lagIdx = 1:length(Params.FuncConLagval)
    lagVal = Params.FuncConLagval(lagIdx);
    lagFieldName = sprintf('adjM%.fmslag', lagVal);
    
    if ~isfield(expData.NetMet, lagFieldName)
        continue;
    end
    
    cartographyData = struct();
    
    % Node cartography metrics
    cartographyFields = {
        'Z',                % Within-module degree z-score
        'PC',               % Participation coefficient  
        'NCpn1',            % Peripheral nodes
        'NCpn2',            % Non-hub connectors
        'NCpn3',            % Non-hub kinless nodes
        'NCpn4',            % Provincial hubs
        'NCpn5',            % Connector hubs
        'NCpn6',            % Kinless hubs
        'Ci'                % Community affiliation
    };
    
    netMetLag = expData.NetMet.(lagFieldName);
    
    for i = 1:length(cartographyFields)
        fieldName = cartographyFields{i};
        if isfield(netMetLag, fieldName)
            cartographyData.(fieldName) = netMetLag.(fieldName);
        end
    end
    
    % Add boundary parameters (needed for node classification)
    cartographyData.boundaries = struct();
    boundaryNames = {'hubBoundaryWMdDeg', 'periPartCoef', 'proHubpartCoef', ...
                     'nonHubconnectorPartCoef', 'connectorHubPartCoef'};
    
    % Check for lag-specific boundaries first
    for i = 1:length(boundaryNames)
        boundaryField = sprintf('%s_%.fmsLag', boundaryNames{i}, lagVal);
        if isfield(Params, boundaryField)
            cartographyData.boundaries.(boundaryNames{i}) = Params.(boundaryField);
        elseif isfield(Params, boundaryNames{i})
            cartographyData.boundaries.(boundaryNames{i}) = Params.(boundaryNames{i});
        end
    end
    
    % Add metadata
    cartographyData.lagValue = lagVal;
    cartographyData.expName = expName;
    cartographyData.Grp = expData.Info.Grp;
    cartographyData.DIV = expData.Info.DIV;
    
    % Add coordinates and channels
    if isfield(expData, 'coords')
        cartographyData.coords = expData.coords;
    end
    if isfield(expData, 'channels')
        cartographyData.channels = expData.channels;
    elseif isfield(expData.Info, 'channels')
        cartographyData.channels = expData.Info.channels;
    end
    
    % Save cartography data
    filename = sprintf('%s_nodeCartography_lag%.f.mat', expName, lagVal);
    filePath = fullfile(saveFolder, filename);
    save(filePath, 'cartographyData');
    fprintf('Exported node cartography data (lag %.f): %s\n', lagVal, filePath);
end
end