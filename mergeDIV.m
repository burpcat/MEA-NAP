function mergeDIV(divFolders, outputFolder)
% Syntax:
%   mergeDIV(divFolders, outputFolder)
%
% Inputs:
%   divFolders - Cell array of strings containing paths to input DIV folders
%   outputFolder - String containing the path to the output folder
%
% Example:
%   mergeDIV({'path/to/DIV38', 'path/to/DIV55'}, 'path/to/merger3855_nap')
%

    % Display start message
    fprintf('Starting DIV folder merging process...\n');
    
    % Step 1: Create the folder structure
    createFolderStructure(divFolders, outputFolder);
    
    % Step 2: Copy files that remain the same
    copyUnchangedFiles(divFolders, outputFolder);
    
    % Step 3: Copy data files from all input folders
    copyDataFiles(divFolders, outputFolder);
    
    % Step 4: Merge CSV files
    mergeCSVFiles(divFolders, outputFolder);
    
    % Step 5: Generate Parameter files (CSV and MAT)
    generateParameterFiles(divFolders, outputFolder);
    
    fprintf('DIV merging completed successfully.\n');
end

function createFolderStructure(divFolders, outputFolder)
    % Create the main output folder if it doesn't exist
    if ~exist(outputFolder, 'dir')
        fprintf('Creating output folder: %s\n', outputFolder);
        mkdir(outputFolder);
    else
        fprintf('Output folder already exists: %s\n', outputFolder);
    end
    
    % Create the necessary subfolders based on the first DIV folder structure
    % (assuming all DIV folders have a similar structure)
    if ~isempty(divFolders)
        firstDivFolder = divFolders{1};
        
        % Create 1_SpikeDetection folder
        spikeDetectionFolder = fullfile(outputFolder, '1_SpikeDetection');
        if ~exist(spikeDetectionFolder, 'dir')
            fprintf('Creating folder: %s\n', spikeDetectionFolder);
            mkdir(spikeDetectionFolder);
        end
        
        % Create 1_SpikeDetection subfolders
        spikeDetectedDataFolder = fullfile(spikeDetectionFolder, '1A_SpikeDetectedData');
        if ~exist(spikeDetectedDataFolder, 'dir')
            fprintf('Creating folder: %s\n', spikeDetectedDataFolder);
            mkdir(spikeDetectedDataFolder);
        end
        
        spikeDetectionChecksFolder = fullfile(spikeDetectionFolder, '1B_SpikeDetectionChecks');
        if ~exist(spikeDetectionChecksFolder, 'dir')
            fprintf('Creating folder: %s\n', spikeDetectionChecksFolder);
            mkdir(spikeDetectionChecksFolder);
        end
        
        % Create folders for different genotypes in SpikeDetectionChecks
        genotypeFolders = {'MUT', 'WM5050', 'WM5050A', 'WM8020A', 'WM9505A', 'WT'};
        for i = 1:length(genotypeFolders)
            genotypeFolder = fullfile(spikeDetectionChecksFolder, genotypeFolders{i});
            if ~exist(genotypeFolder, 'dir')
                fprintf('Creating folder: %s\n', genotypeFolder);
                mkdir(genotypeFolder);
            end
        end
        
        % Create ExperimentMatFiles folder
        experimentMatFilesFolder = fullfile(outputFolder, 'ExperimentMatFiles');
        if ~exist(experimentMatFilesFolder, 'dir')
            fprintf('Creating folder: %s\n', experimentMatFilesFolder);
            mkdir(experimentMatFilesFolder);
        end
    else
        warning('No DIV folders provided. Cannot create folder structure.');
    end
end

function copyUnchangedFiles(divFolders, outputFolder)
    % Copy channel_layout.png from the first DIV folder (assuming it's the same in all folders)
    if ~isempty(divFolders)
        firstDivFolder = divFolders{1};
        
        % Copy channel_layout.png
        sourceFile = fullfile(firstDivFolder, 'channel_layout.png');
        destFile = fullfile(outputFolder, 'channel_layout.png');
        
        if exist(sourceFile, 'file')
            fprintf('Copying file: %s to %s\n', sourceFile, destFile);
            copyfile(sourceFile, destFile);
        else
            warning('channel_layout.png not found in %s', firstDivFolder);
        end
    else
        warning('No DIV folders provided. Cannot copy files.');
    end
end

function copyDataFiles(divFolders, outputFolder)
    % Get the output folder base name (used for renaming files)
    [~, outputFolderName] = fileparts(outputFolder);
    
    % Process each input DIV folder
    for i = 1:length(divFolders)
        divFolder = divFolders{i};
        [~, divName] = fileparts(divFolder);
        fprintf('Processing %s...\n', divName);
        
        % 1. Copy files from 1A_SpikeDetectedData
        sourceSpikeDataFolder = fullfile(divFolder, '1_SpikeDetection', '1A_SpikeDetectedData');
        destSpikeDataFolder = fullfile(outputFolder, '1_SpikeDetection', '1A_SpikeDetectedData');
        
        if exist(sourceSpikeDataFolder, 'dir')
            % Get all .mat files
            spikeFiles = dir(fullfile(sourceSpikeDataFolder, '*.mat'));
            
            % Copy each file
            for j = 1:length(spikeFiles)
                sourceFile = fullfile(sourceSpikeDataFolder, spikeFiles(j).name);
                destFile = fullfile(destSpikeDataFolder, spikeFiles(j).name);
                
                if ~exist(destFile, 'file')
                    fprintf('Copying spike data file: %s\n', spikeFiles(j).name);
                    copyfile(sourceFile, destFile);
                else
                    fprintf('Spike data file already exists, skipping: %s\n', spikeFiles(j).name);
                end
            end
        else
            warning('SpikeDetectedData folder not found in %s', divFolder);
        end
        
        % 2. Copy files from 1B_SpikeDetectionChecks
        sourceSpikeChecksFolder = fullfile(divFolder, '1_SpikeDetection', '1B_SpikeDetectionChecks');
        destSpikeChecksFolder = fullfile(outputFolder, '1_SpikeDetection', '1B_SpikeDetectionChecks');
        
        if exist(sourceSpikeChecksFolder, 'dir')
            % Get all genotype folders
            genotypeFolders = dir(sourceSpikeChecksFolder);
            genotypeFolders = genotypeFolders([genotypeFolders.isdir]);  % Keep only directories
            genotypeFolders = genotypeFolders(~ismember({genotypeFolders.name}, {'.', '..'}));  % Remove . and ..
            
            % Process each genotype folder
            for j = 1:length(genotypeFolders)
                genotypeFolder = genotypeFolders(j).name;
                sourceGenotypeFolder = fullfile(sourceSpikeChecksFolder, genotypeFolder);
                destGenotypeFolder = fullfile(destSpikeChecksFolder, genotypeFolder);
                
                % Create destination genotype folder if it doesn't exist
                if ~exist(destGenotypeFolder, 'dir')
                    mkdir(destGenotypeFolder);
                end
                
                % Get all experiment folders within this genotype
                experimentFolders = dir(sourceGenotypeFolder);
                experimentFolders = experimentFolders([experimentFolders.isdir]);  % Keep only directories
                experimentFolders = experimentFolders(~ismember({experimentFolders.name}, {'.', '..'}));  % Remove . and ..
                
                % Process each experiment folder
                for k = 1:length(experimentFolders)
                    experimentFolder = experimentFolders(k).name;
                    sourceExperimentFolder = fullfile(sourceGenotypeFolder, experimentFolder);
                    destExperimentFolder = fullfile(destGenotypeFolder, experimentFolder);
                    
                    % Create destination experiment folder if it doesn't exist
                    if ~exist(destExperimentFolder, 'dir')
                        fprintf('Creating experiment folder: %s\n', destExperimentFolder);
                        mkdir(destExperimentFolder);
                    end
                    
                    % Copy all .png files
                    pngFiles = dir(fullfile(sourceExperimentFolder, '*.png'));
                    for l = 1:length(pngFiles)
                        sourceFile = fullfile(sourceExperimentFolder, pngFiles(l).name);
                        destFile = fullfile(destExperimentFolder, pngFiles(l).name);
                        
                        if ~exist(destFile, 'file')
                            fprintf('Copying PNG file: %s to %s\n', pngFiles(l).name, destExperimentFolder);
                            copyfile(sourceFile, destFile);
                        else
                            fprintf('PNG file already exists, skipping: %s\n', pngFiles(l).name);
                        end
                    end
                end
            end
        else
            warning('SpikeDetectionChecks folder not found in %s', divFolder);
        end
        
        % 3. Copy files from ExperimentMatFiles
        sourceExpMatFolder = fullfile(divFolder, 'ExperimentMatFiles');
        destExpMatFolder = fullfile(outputFolder, 'ExperimentMatFiles');
        
        if exist(sourceExpMatFolder, 'dir')
            % Get all .mat files
            matFiles = dir(fullfile(sourceExpMatFolder, '*.mat'));
            
            % Copy each file and rename with the output folder name
            for j = 1:length(matFiles)
                sourceFile = fullfile(sourceExpMatFolder, matFiles(j).name);
                
                % Replace the DIV name suffix with the output folder name
                [~, baseFileName, ~] = fileparts(matFiles(j).name);
                baseFileName = regexprep(baseFileName, ['_' divName '$'], ['_' outputFolderName]);
                destFile = fullfile(destExpMatFolder, [baseFileName '.mat']);
                
                if ~exist(destFile, 'file')
                    fprintf('Copying and renaming mat file: %s to %s\n', matFiles(j).name, [baseFileName '.mat']);
                    copyfile(sourceFile, destFile);
                else
                    fprintf('Mat file already exists, skipping: %s\n', [baseFileName '.mat']);
                end
            end
        else
            warning('ExperimentMatFiles folder not found in %s', divFolder);
        end
    end
end

function mergeCSVFiles(divFolders, outputFolder)
    % Initialize empty table for merging
    mergedData = [];
    
    % Process each input DIV folder
    for i = 1:length(divFolders)
        divFolder = divFolders{i};
        fprintf('Processing CSV files in %s...\n', divFolder);
        
        % Find all CSV files in the input folder
        csvFiles = dir(fullfile(divFolder, '*.csv'));
        
        % Filter out Parameters.csv files
        nonParameterCSVs = {};
        for j = 1:length(csvFiles)
            if ~contains(csvFiles(j).name, 'Parameters')
                nonParameterCSVs{end+1} = csvFiles(j).name;
            end
        end
        
        % Check if we found exactly one non-Parameters CSV
        if length(nonParameterCSVs) == 1
            csvFile = nonParameterCSVs{1};
            fprintf('Found CSV file to merge: %s\n', csvFile);
            
            % Read the CSV file
            sourceFile = fullfile(divFolder, csvFile);
            try
                % Try reading with readtable (better for CSV with headers)
                csvData = readtable(sourceFile, 'TextType', 'string');
                
                % If this is the first data we're reading, use it to initialize the merged data
                if isempty(mergedData)
                    mergedData = csvData;
                else
                    % Otherwise, concatenate with existing data
                    % Make sure column names match (case-sensitive)
                    if isequal(csvData.Properties.VariableNames, mergedData.Properties.VariableNames)
                        mergedData = [mergedData; csvData];
                    else
                        warning('CSV file %s has different column names than previous files. Skipping.', sourceFile);
                        fprintf('Expected columns: %s\n', strjoin(mergedData.Properties.VariableNames, ', '));
                        fprintf('Found columns: %s\n', strjoin(csvData.Properties.VariableNames, ', '));
                    end
                end
            catch e
                warning('Error reading CSV file %s: %s', sourceFile, e.message);
                
                % Try alternative CSV reading methods if readtable fails
                try
                    % Try csvread (simpler, numeric only, no headers)
                    csvData = csvread(sourceFile);
                    
                    % If this is the first data we're reading, use it
                    if isempty(mergedData)
                        mergedData = csvData;
                    else
                        % Check if dimensions match
                        if size(csvData, 2) == size(mergedData, 2)
                            mergedData = [mergedData; csvData];
                        else
                            warning('CSV file %s has different number of columns than previous files. Skipping.', sourceFile);
                        end
                    end
                catch e2
                    warning('Failed to read CSV file %s using alternative method: %s', sourceFile, e2.message);
                end
            end
        elseif isempty(nonParameterCSVs)
            warning('No non-Parameters CSV files found in %s', divFolder);
        else
            warning('Multiple non-Parameters CSV files found in %s, expected only one.', divFolder);
            fprintf('Found files: %s\n', strjoin(nonParameterCSVs, ', '));
        end
    end
    
    % Write the merged data to the output folder
    if ~isempty(mergedData)
        destFile = fullfile(outputFolder, 'div_merged.csv');
        fprintf('Writing merged CSV data to %s\n', destFile);
        
        try
            if istable(mergedData)
                % Write table to CSV
                writetable(mergedData, destFile);
            else
                % Write matrix to CSV
                csvwrite(destFile, mergedData);
            end
            fprintf('Successfully wrote merged CSV file\n');
        catch e
            warning('Error writing merged CSV file: %s', e.message);
        end
    else
        warning('No CSV data was merged, cannot write output file');
    end
end

function generateParameterFiles(divFolders, outputFolder)
    % Get output folder name
    [~, outputFolderName] = fileparts(outputFolder);
    
    % Initialize arrays to collect parameter data from all input folders
    allParameterCsvData = {};
    allParameterMatData = {};
    parameterCsvFiles = {};
    
    % Process each input DIV folder
    fprintf('Collecting parameter data from input folders...\n');
    for i = 1:length(divFolders)
        divFolder = divFolders{i};
        [~, divName] = fileparts(divFolder);
        
        % Find the Parameters CSV file
        parameterCsvFile = fullfile(divFolder, ['Parameters_' divName '.csv']);
        if exist(parameterCsvFile, 'file')
            fprintf('Found Parameter CSV file: %s\n', parameterCsvFile);
            parameterCsvFiles{end+1} = parameterCsvFile;
            try
                % Read the CSV file
                paramData = readtable(parameterCsvFile, 'TextType', 'string');
                allParameterCsvData{end+1} = paramData;
            catch e
                warning('Error reading Parameter CSV file %s: %s', parameterCsvFile, e.message);
            end
        else
            warning('Parameter CSV file not found in %s', divFolder);
        end
        
        % Find the Parameters MAT file
        parameterMatFile = fullfile(divFolder, ['Parameters_' divName '.mat']);
        if exist(parameterMatFile, 'file')
            fprintf('Found Parameter MAT file: %s\n', parameterMatFile);
            try
                % Load the MAT file
                matData = load(parameterMatFile);
                allParameterMatData{end+1} = matData;
            catch e
                warning('Error loading Parameter MAT file %s: %s', parameterMatFile, e.message);
            end
        else
            warning('Parameter MAT file not found in %s', divFolder);
        end
    end
    
    % Generate merged Parameter CSV file
    fprintf('Generating merged Parameter CSV file...\n');
    if ~isempty(parameterCsvFiles)
        % Create the path to the merged DIV CSV
        divMergedCsvPath = fullfile(outputFolder, 'div_merged.csv');
        
        % Call the function to merge parameter CSV files
        mergeParameterCSV(parameterCsvFiles, outputFolder, divMergedCsvPath);
    else
        warning('No Parameter CSV data found in input folders');
    end
    
    % Generate merged Parameter MAT file
    fprintf('Generating merged Parameter MAT file...\n');
    if ~isempty(allParameterMatData)
        % Output file path
        outputMatFile = fullfile(outputFolder, ['Parameters_' outputFolderName '.mat']);
        fprintf('Will write the parameter MAT file to: %s\n', outputMatFile);
        
        % TODO: Implement the MAT file generation logic
        % This section is intentionally left empty as requested
        fprintf('MAT file generation logic not yet implemented\n');
    else
        warning('No Parameter MAT data found in input folders');
    end
end

function mergeParameterCSV(parameterCsvFiles, outputFolder, divMergedCsvPath)
    % mergeParameterCSV - Merges parameter CSV files from multiple DIV folders
    %
    % Syntax:
    %   mergeParameterCSV(parameterCsvFiles, outputFolder, divMergedCsvPath)
    %
    % Inputs:
    %   parameterCsvFiles - Cell array of strings containing paths to parameter CSV files
    %   outputFolder - String containing the path to the output folder
    %   divMergedCsvPath - Path to the merged DIV CSV file
    
    % Get output folder name
    [~, outputFolderName] = fileparts(outputFolder);
    
    % Check if we have any input files
    if isempty(parameterCsvFiles)
        error('No parameter CSV files provided');
    end
    
    % Initialize channel and coordinate offsets
    channelOffset = 0;
    coordOffset = 0;
    
    % Process the first file as our base
    fprintf('Processing first parameter file as base: %s\n', parameterCsvFiles{1});
    mergedData = readtable(parameterCsvFiles{1}, 'TextType', 'string');
    
    % Update the specific fields in the merged data
    if ismember('outputDataFolderName', mergedData.Properties.VariableNames)
        mergedData.outputDataFolderName(1) = string(outputFolderName);
        fprintf('Updated outputDataFolderName to: %s\n', outputFolderName);
    end
    
    if ismember('rawData', mergedData.Properties.VariableNames)
        mergedData.rawData(1) = "default_value";
        fprintf('Set rawData to default value\n');
    end
    
    if ismember('spreadSheetFileName', mergedData.Properties.VariableNames)
        mergedData.spreadSheetFileName(1) = string(divMergedCsvPath);
        fprintf('Updated spreadSheetFileName to: %s\n', divMergedCsvPath);
    end
    
    % Find the highest channel and coordinate indices in the first file
    channelColPattern = 'channels_(\d+)_(\d+)';
    coordColPattern = 'coords_(\d+)_(\d+)';
    varNames = mergedData.Properties.VariableNames;
    
    maxChannelRow = 0;
    maxCoordRow = 0;
    
    for i = 1:length(varNames)
        % Check if column is a channel column
        channelMatch = regexp(varNames{i}, channelColPattern, 'tokens');
        if ~isempty(channelMatch)
            row = str2double(channelMatch{1}{1});
            maxChannelRow = max(maxChannelRow, row);
        end
        
        % Check if column is a coord column
        coordMatch = regexp(varNames{i}, coordColPattern, 'tokens');
        if ~isempty(coordMatch)
            row = str2double(coordMatch{1}{1});
            maxCoordRow = max(maxCoordRow, row);
        end
    end
    
    fprintf('First file has channels up to row %d and coords up to row %d\n', maxChannelRow, maxCoordRow);
    
    % Update offsets for next file
    channelOffset = maxChannelRow;
    coordOffset = maxCoordRow;
    
    % If we have a DivNm column, we'll need to rename it if there are multiple files
    hasDivNm = ismember('DivNm', mergedData.Properties.VariableNames);
    
    % Process each additional file
    for fileIdx = 2:length(parameterCsvFiles)
        fprintf('Processing parameter file %d of %d: %s\n', fileIdx, length(parameterCsvFiles), parameterCsvFiles{fileIdx});
        
        % Read the current file
        currentFile = readtable(parameterCsvFiles{fileIdx}, 'TextType', 'string');
        currentVarNames = currentFile.Properties.VariableNames;
        
        % Find channel and coordinate columns
        channelCols = {};
        coordCols = {};
        
        for i = 1:length(currentVarNames)
            varName = currentVarNames{i};
            
            % Check if column is a channel column
            channelMatch = regexp(varName, channelColPattern, 'tokens');
            if ~isempty(channelMatch)
                channelCols{end+1} = varName;
            end
            
            % Check if column is a coord column
            coordMatch = regexp(varName, coordColPattern, 'tokens');
            if ~isempty(coordMatch)
                coordCols{end+1} = varName;
            end
        end
        
        % Add channel columns with renamed indices to merged data
        for i = 1:length(channelCols)
            oldName = channelCols{i};
            channelMatch = regexp(oldName, channelColPattern, 'tokens');
            row = str2double(channelMatch{1}{1});
            col = str2double(channelMatch{1}{2});
            
            % Create new column name with adjusted row index
            newRow = row + channelOffset;
            newName = sprintf('channels_%d_%d', newRow, col);
            
            % Add this column to mergedData
            mergedData.(newName) = currentFile.(oldName)(1);
        end
        
        % Add coordinate columns with renamed indices to merged data
        for i = 1:length(coordCols)
            oldName = coordCols{i};
            coordMatch = regexp(oldName, coordColPattern, 'tokens');
            row = str2double(coordMatch{1}{1});
            col = str2double(coordMatch{1}{2});
            
            % Create new column name with adjusted row index
            newRow = row + coordOffset;
            newName = sprintf('coords_%d_%d', newRow, col);
            
            % Add this column to mergedData
            mergedData.(newName) = currentFile.(oldName)(1);
        end
        
        % Handle DivNm column
        if fileIdx == 2 && hasDivNm
            % For the second file, rename DivNm to DivNm1 in the merged data
            divNmValue = mergedData.DivNm(1);
            
            % Create a new column DivNm1 with the value from DivNm
            mergedData.DivNm1 = divNmValue;
            
            % Remove the original DivNm column
            mergedData = removevars(mergedData, 'DivNm');
        end
        
        % Add DivNm from current file as DivNmX
        if ismember('DivNm', currentFile.Properties.VariableNames)
            newDivNmCol = sprintf('DivNm%d', fileIdx);
            mergedData.(newDivNmCol) = currentFile.DivNm(1);
        end
        
        % Update offsets for next file (assuming standard 24 rows per file)
        channelOffset = channelOffset + 24;
        coordOffset = coordOffset + 24;
    end
    
    % Write the final merged data to output file
    outputCsvFile = fullfile(outputFolder, ['Parameters_' outputFolderName '.csv']);
    writetable(mergedData, outputCsvFile);
    fprintf('Successfully wrote merged parameter CSV file to: %s\n', outputCsvFile);
end