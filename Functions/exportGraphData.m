function exportGraphData(dataStage, expData, Params, ExN)
% Centralized function to export graph data at different pipeline stages
%
% Inputs:
%   dataStage - string indicating the pipeline stage ('spike', 'adjmatrix', 'netmet', etc.)
%   expData - structure containing experiment data (can include Info, adjMs, NetMet, etc.)
%   Params - pipeline parameters structure
%   ExN - experiment number (optional, only needed for certain stages)
%
% This function automatically determines what data is available and saves it
% to the GraphData folder for later visualization.

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

% Create folder structure if needed
groupFolder = fullfile(Params.graphDataFolder, groupName);
if ~isfolder(groupFolder)
    mkdir(groupFolder);
end

expFolder = fullfile(groupFolder, expName);
if ~isfolder(expFolder)
    mkdir(expFolder);
end

% Select what to export based on the data stage
switch lower(dataStage)
    case 'spike'
        % Export spike detection data
        exportSpikeData(expData, expName, expFolder);
        
    case 'ephys'
        % Export neuronal activity data (Ephys)
        exportEphysData(expData, expName, expFolder, Params, ExN);
        
    case 'adjmatrix'
        % Export adjacency matrices
        exportAdjMatrixData(expData, expName, expFolder, Params);
        
    case 'netmet'
        % Export network metrics
        exportNetMetData(expData, expName, expFolder, Params);
        
    case 'nodecartography'
        % Export node cartography data
        exportNodeCartographyData(expData, expName, expFolder, Params);
        
    otherwise
        warning('Unknown data stage: %s', dataStage);
end

end

function exportSpikeData(expData, expName, saveFolder)
% Export spike detection data

% Check if we have spike data
if isfield(expData, 'spikeTimes')
    % Prepare spike data
    spikeData = struct();
    spikeData.spikeTimes = expData.spikeTimes;
    
    if isfield(expData, 'channels')
        spikeData.channels = expData.channels;
    elseif isfield(expData, 'Info') && isfield(expData.Info, 'channels')
        spikeData.channels = expData.Info.channels;
    end
    
    % Save spike data
    filePath = fullfile(saveFolder, [expName, '_spikeData.mat']);
    save(filePath, 'spikeData');
    fprintf('Saved spike data: %s\n', filePath);
end

% Check if we have spike matrix
if isfield(expData, 'spikeMatrix')
    % Prepare raster data
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
    
    % Save raster data
    filePath = fullfile(saveFolder, [expName, '_spikeRaster.mat']);
    save(filePath, 'rasterData');
    fprintf('Saved spike raster data: %s\n', filePath);
end
end

function exportEphysData(expData, expName, saveFolder, Params, ExN)
Export neuronal activity data (Ephys)

if ~isfield(expData, 'Ephys')
    return;
end

Prepare activity data
activityData = struct();
ephysFields = {'FR', 'channelBurstRate', 'channelBurstDur', 'channelFracSpikesInBursts', ...
    'channelISIwithinBurst', 'channeISIoutsideBurst'};

for i = 1:length(ephysFields)
    if isfield(expData.Ephys, ephysFields{i})
        activityData.(ephysFields{i}) = expData.Ephys.(ephysFields{i});
    end
end

Add coordinates and channels
if isfield(expData, 'coords')
    activityData.coords = expData.coords;
elseif isfield(Params, 'coords') && ExN <= length(Params.coords)
    activityData.coords = Params.coords{ExN};
end

if isfield(expData, 'channels')
    activityData.channels = expData.channels;
elseif isfield(expData.Info, 'channels')
    activityData.channels = expData.Info.channels;
end

Save activity data
filePath = fullfile(saveFolder, [expName, '_electrodeSpikeActivity.mat']);
save(filePath, 'activityData');
fprintf('Saved electrode activity data: %s\n', filePath);
end

function exportAdjMatrixData(expData, expName, saveFolder, Params)
% Export adjacency matrices

if ~isfield(expData, 'adjMs')
    return;
end

% Export adjacency matrices for each lag value
for lagIdx = 1:length(Params.FuncConLagval)
    lagVal = Params.FuncConLagval(lagIdx);
    adjMatrixFieldName = sprintf('adjM%.fmslag', lagVal);
    
    if isfield(expData.adjMs, adjMatrixFieldName)
        adjM = expData.adjMs.(adjMatrixFieldName);
        
        % Save adjacency matrix
        filename = sprintf('%s_adjMatrix_lag%d.mat', expName, lagVal);
        filePath = fullfile(saveFolder, filename);
        save(filePath, 'adjM');
        fprintf('Saved adjacency matrix (lag %d): %s\n', lagVal, filePath);
    end
end
end

function exportNetMetData(expData, expName, saveFolder, Params)
% Export network metrics

if ~isfield(expData, 'NetMet')
    return;
end

% Export network metrics for each lag value
for lagIdx = 1:length(Params.FuncConLagval)
    lagVal = Params.FuncConLagval(lagIdx);
    
    % Prepare node-level metrics
    nodeMetrics = struct();
    metricFields = {'degree', 'strength', 'clustering', 'betweenness', 'efficiency_local', ...
        'participation', 'control_average', 'control_modal'};
    
    for mIdx = 1:length(metricFields)
        fieldName = metricFields{mIdx};
        metricFieldName = sprintf('%s%.fmslag', fieldName, lagVal);
        if isfield(expData.NetMet, metricFieldName)
            nodeMetrics.(fieldName) = expData.NetMet.(metricFieldName);
        end
    end
    
    % Add node coordinates and channels
    if isfield(expData, 'coords')
        nodeMetrics.coords = expData.coords;
    end
    if isfield(expData, 'channels')
        nodeMetrics.channels = expData.channels;
    elseif isfield(expData.Info, 'channels')
        nodeMetrics.channels = expData.Info.channels;
    end
    
    % Save node metrics
    nodeFilename = sprintf('%s_nodeMetrics_lag%d.mat', expName, lagVal);
    nodeFilePath = fullfile(saveFolder, nodeFilename);
    save(nodeFilePath, 'nodeMetrics');
    fprintf('Saved node metrics (lag %d): %s\n', lagVal, nodeFilePath);
    
    % Prepare network-level metrics
    netMetrics = struct();
    netMetricFields = {'density', 'efficiency_global', 'modularity', 'smallworldness'};
    
    for mIdx = 1:length(netMetricFields)
        fieldName = netMetricFields{mIdx};
        metricFieldName = sprintf('%s%.fmslag', fieldName, lagVal);
        if isfield(expData.NetMet, metricFieldName)
            netMetrics.(fieldName) = expData.NetMet.(metricFieldName);
        end
    end
    
    % Save network metrics
    netFilename = sprintf('%s_networkMetrics_lag%d.mat', expName, lagVal);
    netFilePath = fullfile(saveFolder, netFilename);
    save(netFilePath, 'netMetrics');
    fprintf('Saved network metrics (lag %d): %s\n', lagVal, netFilePath);
end
end

function exportNodeCartographyData(expData, expName, saveFolder, Params)
% Export node cartography data

if ~isfield(expData, 'NetMet')
    return;
end

% Export node cartography for each lag value
for lagIdx = 1:length(Params.FuncConLagval)
    lagVal = Params.FuncConLagval(lagIdx);
    
    % Check if we have node cartography data for this lag
    hasCartography = false;
    cartographyFields = {'z', 'p', 'nodal_roles'};
    
    for cIdx = 1:length(cartographyFields)
        fieldName = cartographyFields{cIdx};
        metricFieldName = sprintf('%s%.fmslag', fieldName, lagVal);
        if isfield(expData.NetMet, metricFieldName)
            hasCartography = true;
            break;
        end
    end
    
    if ~hasCartography
        continue;
    end
    
    % Prepare cartography data
    cartographyData = struct();
    
    for cIdx = 1:length(cartographyFields)
        fieldName = cartographyFields{cIdx};
        metricFieldName = sprintf('%s%.fmslag', fieldName, lagVal);
        if isfield(expData.NetMet, metricFieldName)
            cartographyData.(fieldName) = expData.NetMet.(metricFieldName);
        end
    end
    
    % Add nodes and roles by region
    roleFieldName = sprintf('NCpn%.fmslag', lagVal);
    if isfield(expData.NetMet, roleFieldName)
        cartographyData.roles = expData.NetMet.(roleFieldName);
    end
    
    % Add boundary values
    cartographyData.boundaries = struct();
    boundaryNames = {'hubBoundaryWMdDeg', 'periPartCoef', 'proHubpartCoef', ...
        'nonHubconnectorPartCoef', 'connectorHubPartCoef'};
    
    for bIdx = 1:length(boundaryNames)
        boundaryFieldName = strcat(boundaryNames{bIdx}, sprintf('_%.fmsLag', lagVal));
        if isfield(Params, boundaryFieldName)
            cartographyData.boundaries.(boundaryNames{bIdx}) = Params.(boundaryFieldName);
        end
    end
    
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
    filename = sprintf('%s_nodeCartography_lag%d.mat', expName, lagVal);
    filePath = fullfile(saveFolder, filename);
    save(filePath, 'cartographyData');
    fprintf('Saved node cartography data (lag %d): %s\n', lagVal, filePath);
end
end