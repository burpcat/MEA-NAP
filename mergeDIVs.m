%% MEA-NAP DIV Merger - Comprehensive DIV Data Merging Script
% This script provides a robust solution for merging multiple DIV (Days In Vitro)
% datasets processed by MEA-NAP. It handles various edge cases, validates data
% compatibility, provides detailed progress feedback, and creates standardized outputs.
%
% Author: Avinash
% Version: 1.0
% Date: 2025-03-22

%% Configuration - EDIT THESE PARAMETERS
% List all DIV folders to merge (path to directories containing the processed DIV data)
DIVfolders = {
    '/Users/avinash/Documents/labwork/MEA-NAP/output/DIV5/div55',
    % Add more DIV folders here
    % '/Users/avinash/Documents/labwork/MEA-NAP/output/DIV10/div10',
    % '/Users/avinash/Documents/labwork/MEA-NAP/output/DIV15/div15'
};

% Output directory for merged results
outputFolder = '/Users/avinash/Documents/labwork/MEA-NAP/output/CombinedDIVs';

% Merge options
options = struct(...
    'createPlots', true, ...                % Generate plots for combined data
    'saveMergedMatFiles', true, ...         % Save merged MAT files
    'exportCSV', true, ...                  % Export data to CSV files
    'overwrite', true, ...                  % Overwrite existing output
    'verbose', true, ...                    % Detailed progress output
    'harmonizeParams', true, ...            % Try to harmonize different parameter structures
    'handleMissingMetrics', 'ignore', ...   % 'ignore', 'error', or 'fillna'
    'parallelProcessing', true, ...         % Use parallel processing if available
    'metricSelectionStrategy', 'union' ...  % 'union', 'intersection', or 'specified'
);

% Specify custom metrics (only used if metricSelectionStrategy = 'specified')
customNetworkMetrics = {}; % Leave empty to use default metrics
customNodeMetrics = {};    % Leave empty to use default metrics

%% Setup and Validation
fprintf('MEA-NAP DIV Merger - Starting\n');
fprintf('=============================\n\n');

% Add MEA-NAP to path if needed
try
    % Try to locate key MEA-NAP functions to check if they're on the path
    if ~exist('combineExpNetworkData', 'file') || ~exist('compileAllExpData', 'file')
        % Try to find MEA-NAP directory
        currentFile = mfilename('fullpath');
        [currentDir, ~, ~] = fileparts(currentFile);
        if exist(fullfile(currentDir, 'MEApipeline.m'), 'file')
            % We're in the MEA-NAP directory
            addpath(genpath(currentDir));
            fprintf('✓ Added MEA-NAP directory to path\n');
        else
            % Try parent directory
            [parentDir, ~, ~] = fileparts(currentDir);
            if exist(fullfile(parentDir, 'MEApipeline.m'), 'file')
                addpath(genpath(parentDir));
                fprintf('✓ Added MEA-NAP directory to path\n');
            else
                warning('Could not find MEA-NAP directory automatically. Please ensure MEA-NAP is on your MATLAB path.');
            end
        end
    end
catch
    warning('Error adding MEA-NAP to path. Ensure MEA-NAP functions are accessible.');
end

% Check if DIV folders exist
validFolders = {};
for i = 1:length(DIVfolders)
    if ~exist(DIVfolders{i}, 'dir')
        warning('DIV folder does not exist: %s. Skipping.', DIVfolders{i});
    else
        % Check for essential files
        paramFiles = dir(fullfile(DIVfolders{i}, 'Parameters_*.mat'));
        if isempty(paramFiles)
            warning('No parameter file found in folder: %s. Skipping.', DIVfolders{i});
        else
            if ~exist(fullfile(DIVfolders{i}, 'ExperimentMatFiles'), 'dir')
                warning('ExperimentMatFiles folder missing in: %s. Skipping.', DIVfolders{i});
            else
                validFolders{end+1} = DIVfolders{i};
                fprintf('✓ Validated DIV folder: %s\n', DIVfolders{i});
            end
        end
    end
end

if isempty(validFolders)
    error('No valid DIV folders found. Please check your folder paths.');
end

% Create output directory if it doesn't exist
if ~exist(outputFolder, 'dir')
    if options.verbose, fprintf('Creating output directory: %s\n', outputFolder); end
    [success, msg] = mkdir(outputFolder);
    if ~success
        error('Failed to create output folder: %s. Error: %s', outputFolder, msg);
    end
else
    if options.verbose, fprintf('Output directory already exists: %s\n', outputFolder); end
    if ~options.overwrite
        error('Output directory already exists and overwrite is set to false.');
    end
end

%% Load parameter files and validate compatibility
if options.verbose, fprintf('\nLoading and validating parameter files...\n'); end

% Load all parameter files
paramStructs = cell(1, length(validFolders));
for i = 1:length(validFolders)
    paramFiles = dir(fullfile(validFolders{i}, 'Parameters_*.mat'));
    [~, idx] = max([paramFiles.datenum]); % Get most recent
    paramFile = fullfile(validFolders{i}, paramFiles(idx).name);
    if options.verbose, fprintf('Loading parameters from: %s\n', paramFile); end
    
    try
        tempData = load(paramFile);
        if isfield(tempData, 'Params')
            paramStructs{i} = tempData.Params;
        else
            warning('No Params variable found in %s. Trying to find parameter structure.', paramFile);
            fnames = fieldnames(tempData);
            for j = 1:length(fnames)
                if isstruct(tempData.(fnames{j}))
                    paramStructs{i} = tempData.(fnames{j});
                    break;
                end
            end
            if ~exist('paramStructs{i}', 'var')
                error('Could not find parameter structure in %s', paramFile);
            end
        end
    catch ME
        warning('Error loading %s: %s', paramFile, ME.message);
        paramStructs{i} = [];
    end
end

% Remove any empty parameter structures
paramStructs = paramStructs(~cellfun(@isempty, paramStructs));
if isempty(paramStructs)
    error('Failed to load any parameter files.');
end

% Check compatibility and harmonize parameters if needed
if options.harmonizeParams && length(paramStructs) > 1
    if options.verbose, fprintf('Harmonizing parameters across multiple DIV datasets...\n'); end
    
    % Create a unified parameter structure based on the most complete one
    paramFieldCounts = cellfun(@(x) numel(fieldnames(x)), paramStructs);
    [~, mostCompleteIdx] = max(paramFieldCounts);
    mergedParams = paramStructs{mostCompleteIdx};
    
    % Find common network metrics
    allNetworkMetrics = {};
    allNodeMetrics = {};
    
    for i = 1:length(paramStructs)
        if isfield(paramStructs{i}, 'networkLevelNetMetToPlot') && ~isempty(paramStructs{i}.networkLevelNetMetToPlot)
            allNetworkMetrics = [allNetworkMetrics, paramStructs{i}.networkLevelNetMetToPlot];
        end
        
        if isfield(paramStructs{i}, 'unitLevelNetMetToPlot') && ~isempty(paramStructs{i}.unitLevelNetMetToPlot)
            allNodeMetrics = [allNodeMetrics, paramStructs{i}.unitLevelNetMetToPlot];
        end
    end
    
    % Remove duplicates
    allNetworkMetrics = unique(allNetworkMetrics);
    allNodeMetrics = unique(allNodeMetrics);
    
    % Set metrics based on strategy
    switch options.metricSelectionStrategy
        case 'union'
            % Use all metrics found in any file
            mergedParams.networkLevelNetMetToPlot = allNetworkMetrics;
            mergedParams.unitLevelNetMetToPlot = allNodeMetrics;
            
        case 'intersection'
            % Use only metrics common to all files
            commonNetworkMetrics = allNetworkMetrics;
            commonNodeMetrics = allNodeMetrics;
            
            for i = 1:length(paramStructs)
                if isfield(paramStructs{i}, 'networkLevelNetMetToPlot') && ~isempty(paramStructs{i}.networkLevelNetMetToPlot)
                    commonNetworkMetrics = intersect(commonNetworkMetrics, paramStructs{i}.networkLevelNetMetToPlot);
                end
                
                if isfield(paramStructs{i}, 'unitLevelNetMetToPlot') && ~isempty(paramStructs{i}.unitLevelNetMetToPlot)
                    commonNodeMetrics = intersect(commonNodeMetrics, paramStructs{i}.unitLevelNetMetToPlot);
                end
            end
            
            mergedParams.networkLevelNetMetToPlot = commonNetworkMetrics;
            mergedParams.unitLevelNetMetToPlot = commonNodeMetrics;
            
        case 'specified'
            % Use user-specified metrics
            if ~isempty(customNetworkMetrics)
                mergedParams.networkLevelNetMetToPlot = customNetworkMetrics;
            end
            if ~isempty(customNodeMetrics)
                mergedParams.unitLevelNetMetToPlot = customNodeMetrics;
            end
    end
    
    % Set paths for output
    mergedParams.outputDataFolder = outputFolder;
    mergedParams.outputDataFolderName = 'CombinedDIVs';
    
    % Save parameters for reference
    paramFile = fullfile(outputFolder, 'Parameters_Merged.mat');
    save(paramFile, 'mergedParams');
    if options.verbose, fprintf('✓ Saved merged parameters to: %s\n', paramFile); end
    
    % Use the merged parameters
    Params = mergedParams;
else
    % Just use parameters from the first valid DIV
    Params = paramStructs{1};
    Params.outputDataFolder = outputFolder;
    Params.outputDataFolderName = 'CombinedDIVs';
end

%% Setup for parallel processing
if options.parallelProcessing
    try
        % Check if Parallel Computing Toolbox is available
        if license('test', 'Distrib_Computing_Toolbox')
            if options.verbose, fprintf('Setting up parallel processing...\n'); end
            if isempty(gcp('nocreate'))
                parpool('local');
                if options.verbose, fprintf('✓ Started parallel pool\n'); end
            else
                if options.verbose, fprintf('✓ Using existing parallel pool\n'); end
            end
        else
            if options.verbose, fprintf('Parallel Computing Toolbox not available. Using serial processing.\n'); end
            options.parallelProcessing = false;
        end
    catch
        options.parallelProcessing = false;
        warning('Error initializing parallel pool. Using serial processing.');
    end
end

%% Combine data across DIVs
if options.verbose, fprintf('\nCombining data across DIVs...\n'); end

% Set up network and node metrics to extract
if ~isfield(Params, 'networkLevelNetMetToPlot') || isempty(Params.networkLevelNetMetToPlot)
    NetMetricsE = {'density', 'meanDegree', 'meanEfficiency', 'modularity', 'smallWorldness'};
    warning('Using default network metrics: %s', strjoin(NetMetricsE, ', '));
else
    NetMetricsE = Params.networkLevelNetMetToPlot;
end

if ~isfield(Params, 'unitLevelNetMetToPlot') || isempty(Params.unitLevelNetMetToPlot)
    NetMetricsC = {'nodeDegree', 'betweennessCentrality', 'participationCoef'};
    warning('Using default node metrics: %s', strjoin(NetMetricsC, ', '));
else
    NetMetricsC = Params.unitLevelNetMetToPlot;
end

% Check for ExperimentMatFiles in each DIV folder
matFileFolders = cell(1, length(validFolders));
for i = 1:length(validFolders)
    expMatFolder = fullfile(validFolders{i}, 'ExperimentMatFiles');
    if exist(expMatFolder, 'dir')
        matFileFolders{i} = expMatFolder;
    else
        warning('ExperimentMatFiles folder not found in %s', validFolders{i});
        matFileFolders{i} = '';
    end
end

% Filter out empty entries
validFolders = validFolders(~cellfun(@isempty, matFileFolders));
matFileFolders = matFileFolders(~cellfun(@isempty, matFileFolders));

if isempty(validFolders)
    error('No valid DIV folders with ExperimentMatFiles found.');
end

% Create combined ExperimentMatFiles folder in output directory
combinedMatFolder = fullfile(outputFolder, 'ExperimentMatFiles');
if ~exist(combinedMatFolder, 'dir')
    mkdir(combinedMatFolder);
end

% Try to extract DIV numbers from folder names for ordering
divNumbers = zeros(size(validFolders));
for i = 1:length(validFolders)
    [~, folderName] = fileparts(validFolders{i});
    matches = regexp(folderName, 'DIV(\d+)', 'tokens');
    if ~isempty(matches) && ~isempty(matches{1})
        divNumbers(i) = str2double(matches{1}{1});
    else
        matches = regexp(folderName, 'div(\d+)', 'tokens');
        if ~isempty(matches) && ~isempty(matches{1})
            divNumbers(i) = str2double(matches{1}{1});
        else
            divNumbers(i) = i;  % Default ordering
        end
    end
end

% Sort DIV folders by DIV number if possible
if ~all(divNumbers == 0)
    [~, sortIdx] = sort(divNumbers);
    validFolders = validFolders(sortIdx);
    matFileFolders = matFileFolders(sortIdx);
    divNumbers = divNumbers(sortIdx);
    
    if options.verbose
        fprintf('DIV folders sorted by DIV number:\n');
        for i = 1:length(validFolders)
            fprintf('  DIV %d: %s\n', divNumbers(i), validFolders{i});
        end
    end
end

% Update Params with DIV information if available
if all(divNumbers > 0)
    Params.DivNm = divNumbers;
end

% Combine data using MEA-NAP functions
try
    fprintf('\nCombining network data...\n');
    
    % First try with combineExpNetworkData
    try
        combinedData = combineExpNetworkData(validFolders, Params, NetMetricsE, NetMetricsC, 'ExperimentMatFiles');
        fprintf('✓ Successfully combined data using combineExpNetworkData\n');
    catch ME
        warning('Error with combineExpNetworkData: %s', ME.message);
        fprintf('Attempting alternative approach with compileAllExpData...\n');
        
        % If that fails, try with compileAllExpData
        networkData = compileAllExpData(validFolders, 'ExperimentMatFiles', Params, 'network');
        nodeData = compileAllExpData(validFolders, 'ExperimentMatFiles', Params, 'node');
        
        % Manually create combinedData structure
        combinedData = struct();
        
        % Assuming compileAllExpData returns data in a format similar to allExpData.WT.TP1.density
        if ~isempty(networkData)
            groupNames = fieldnames(networkData);
            for g = 1:length(groupNames)
                groupName = groupNames{g};
                combinedData.(groupName) = networkData.(groupName);
                
                % Add node-level data if available
                if ~isempty(nodeData) && isfield(nodeData, groupName)
                    timePoints = fieldnames(nodeData.(groupName));
                    for t = 1:length(timePoints)
                        tpName = timePoints{t};
                        if isfield(combinedData.(groupName), tpName)
                            nodeMetricNames = fieldnames(nodeData.(groupName).(tpName));
                            for n = 1:length(nodeMetricNames)
                                nmName = nodeMetricNames{n};
                                combinedData.(groupName).(tpName).(nmName) = nodeData.(groupName).(tpName).(nmName);
                            end
                        end
                    end
                end
            end
        else
            error('Failed to compile network data.');
        end
        
        fprintf('✓ Successfully combined data using compileAllExpData\n');
    end
    
    % Save combined data
    if options.saveMergedMatFiles
        combinedDataFile = fullfile(outputFolder, 'CombinedDIVData.mat');
        save(combinedDataFile, 'combinedData', 'Params', '-v7.3');
        fprintf('✓ Saved combined data to %s\n', combinedDataFile);
    end
    
catch ME
    error('Failed to combine DIV data. Error: %s', ME.message);
end

%% Generate plots
if options.createPlots
    fprintf('\nGenerating plots...\n');
    try
        PlotNetMet(combinedData, Params);
        fprintf('✓ Successfully generated plots\n');
    catch ME
        warning('Error generating plots: %s', ME.message);
        
        % Try alternative plotting
        fprintf('Attempting alternative plotting approach...\n');
        try
            % Get field names from combinedData
            groupNames = fieldnames(combinedData);
            
            for g = 1:length(groupNames)
                group = groupNames{g};
                if isstruct(combinedData.(group))
                    timepoints = fieldnames(combinedData.(group));
                    
                    % Plot some key metrics if they exist
                    try
                        % Create figure directory
                        figDir = fullfile(outputFolder, 'SimplePlots');
                        if ~exist(figDir, 'dir')
                            mkdir(figDir);
                        end
                        
                        % Check for common metrics
                        commonMetrics = {'density', 'meanDegree', 'efficiency', 'clusterCoef', 'modularity'};
                        for m = 1:length(commonMetrics)
                            metric = commonMetrics{m};
                            
                            % Check if this metric exists across timepoints
                            metricExists = true;
                            for t = 1:length(timepoints)
                                if ~isfield(combinedData.(group).(timepoints{t}), metric)
                                    metricExists = false;
                                    break;
                                end
                            end
                            
                            if metricExists
                                % Extract data for this metric across timepoints
                                metricData = zeros(1, length(timepoints));
                                for t = 1:length(timepoints)
                                    metricData(t) = mean(combinedData.(group).(timepoints{t}).(metric), 'omitnan');
                                end
                                
                                % Create simple plot
                                fig = figure('Visible', 'off');
                                plot(1:length(timepoints), metricData, 'o-', 'LineWidth', 2);
                                xticks(1:length(timepoints));
                                if all(divNumbers > 0)
                                    xticklabels(arrayfun(@(x) sprintf('DIV %d', x), divNumbers, 'UniformOutput', false));
                                else
                                    xticklabels(timepoints);
                                end
                                title(sprintf('%s - %s', group, strrep(metric, '_', ' ')));
                                ylabel(strrep(metric, '_', ' '));
                                xlabel('Timepoint');
                                grid on;
                                saveas(fig, fullfile(figDir, sprintf('%s_%s.png', group, metric)));
                                close(fig);
                                
                                fprintf('  Created plot for %s - %s\n', group, metric);
                            end
                        end
                    catch ME2
                        warning('Error in alternative plotting: %s', ME2.message);
                    end
                end
            end
            
            fprintf('✓ Created basic plots in %s\n', figDir);
        catch ME3
            warning('Failed alternative plotting too: %s', ME3.message);
        end
    end
end

%% Export data to CSV if requested
if options.exportCSV
    fprintf('\nExporting data to CSV files...\n');
    try
        % Base CSV file paths
        nodeCSV = fullfile(outputFolder, 'Combined_NodeLevel.csv');
        networkCSV = fullfile(outputFolder, 'Combined_NetworkLevel.csv');
        
        % Export network-level data
        networkTable = table();
        nodeTable = table();
        
        groupNames = fieldnames(combinedData);
        for g = 1:length(groupNames)
            group = groupNames{g};
            
            if isstruct(combinedData.(group))
                timepoints = fieldnames(combinedData.(group));
                
                for t = 1:length(timepoints)
                    tp = timepoints{t};
                    if isstruct(combinedData.(group).(tp))
                        % Get DIV number from timepoint name or use position
                        divNum = t;
                        if all(divNumbers > 0)
                            divNum = divNumbers(t);
                        else
                            matches = regexp(tp, 'TP(\d+)', 'tokens');
                            if ~isempty(matches) && ~isempty(matches{1})
                                divNum = str2double(matches{1}{1});
                            end
                        end
                        
                        % Process network metrics
                        metricNames = fieldnames(combinedData.(group).(tp));
                        for m = 1:length(metricNames)
                            metric = metricNames{m};
                            metricData = combinedData.(group).(tp).(metric);
                            
                            % Skip non-numeric or empty data
                            if ~isnumeric(metricData) || isempty(metricData)
                                continue;
                            end
                            
                            % Determine if this is a node-level or network-level metric
                            if size(metricData, 1) > 1 || size(metricData, 2) > 1
                                % Node-level metric (multi-value)
                                [rows, cols] = size(metricData);
                                for i = 1:rows
                                    for j = 1:min(cols, 1)  % Only take first column for simplicity
                                        rowData = table(group, divNum, metric, i, metricData(i,j), ...
                                            'VariableNames', {'Group', 'DIV', 'Metric', 'NodeID', 'Value'});
                                        nodeTable = [nodeTable; rowData];
                                    end
                                end
                            else
                                % Network-level metric (single value)
                                rowData = table(group, divNum, metric, metricData, ...
                                    'VariableNames', {'Group', 'DIV', 'Metric', 'Value'});
                                networkTable = [networkTable; rowData];
                            end
                        end
                    end
                end
            end
        end
        
        % Write tables to CSV
        if ~isempty(networkTable)
            writetable(networkTable, networkCSV);
            fprintf('✓ Saved network-level data to %s\n', networkCSV);
        else
            warning('No network-level data to export.');
        end
        
        if ~isempty(nodeTable)
            writetable(nodeTable, nodeCSV);
            fprintf('✓ Saved node-level data to %s\n', nodeCSV);
        else
            warning('No node-level data to export.');
        end
        
    catch ME
        warning('Error exporting to CSV: %s', ME.message);
    end
end

%% Final report
fprintf('\n=============================\n');
fprintf('MEA-NAP DIV Merger - Complete\n\n');
fprintf('Summary:\n');
fprintf('- Processed %d DIV folders\n', length(validFolders));
fprintf('- Output directory: %s\n', outputFolder);
if options.createPlots
    fprintf('- Generated combined network plots\n');
end
if options.saveMergedMatFiles
    fprintf('- Saved merged MAT files\n');
end
if options.exportCSV
    fprintf('- Exported data to CSV files\n');
end
fprintf('\nUse this data to analyze how network metrics change across DIVs!\n');