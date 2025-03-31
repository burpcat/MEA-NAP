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

    % Read all parameter CSV files
    allData = cell(length(parameterCsvFiles), 1);
    for i = 1:length(parameterCsvFiles)
        % Read table with default handling of empty values
        allData{i} = readtable(parameterCsvFiles{i}, 'TextType', 'string');
    end

    % Create a new table for the merged data, starting with the first file
    mergedData = allData{1};

    % Update the specific fields in the merged data
    if ismember('outputDataFolderName', mergedData.Properties.VariableNames)
        mergedData.outputDataFolderName(1) = string(outputFolderName);
    end

    % Don't modify rawData - keep original value from first CSV
    % if ismember('rawData', mergedData.Properties.VariableNames)
    %     mergedData.rawData(1) = "default_value";
    % end

    if ismember('spreadSheetFileName', mergedData.Properties.VariableNames)
        mergedData.spreadSheetFileName(1) = string(divMergedCsvPath);
    end

    % Define column patterns
    channelColPattern = 'channels_(\d+)_(\d+)';
    coordColPattern = 'coords_(\d+)_(\d+)';

    % 1. First, collect all columns that are not channels or coords
    allNonSpecialCols = {};
    for i = 1:length(allData)
        varNames = allData{i}.Properties.VariableNames;
        for j = 1:length(varNames)
            varName = varNames{j};

            % Skip channel and coord columns
            if ~isempty(regexp(varName, channelColPattern, 'once')) || ...
               ~isempty(regexp(varName, coordColPattern, 'once'))
                continue;
            end
            
            % Skip automatically generated Var columns
            if ~isempty(regexp(varName, '^Var\d+$', 'once'))
                continue;
            end

            % Add non-special column to our list if not already there
            if ~ismember(varName, allNonSpecialCols)
                allNonSpecialCols{end+1} = varName;
            end
        end
    end

    % 2. Create a new table to hold the result
    result = table();

    % 3. Process all channel columns first (in order)
    channelOffset = 0;
    for i = 1:length(allData)
        varNames = allData{i}.Properties.VariableNames;

        % Find all channel columns in this file
        channelCols = {};
        channelRowIndices = [];
        channelColIndices = [];

        for j = 1:length(varNames)
            varName = varNames{j};
            channelMatch = regexp(varName, channelColPattern, 'tokens');

            if ~isempty(channelMatch)
                row = str2double(channelMatch{1}{1});
                col = str2double(channelMatch{1}{2});

                channelCols{end+1} = varName;
                channelRowIndices(end+1) = row;
                channelColIndices(end+1) = col;
            end
        end

        % Sort the channel columns by row, then by column
        [~, sortIdx] = sortrows([channelRowIndices', channelColIndices']);
        channelCols = channelCols(sortIdx);

        % Add channel columns to result, renaming as needed
        for j = 1:length(channelCols)
            oldName = channelCols{j};
            channelMatch = regexp(oldName, channelColPattern, 'tokens');
            row = str2double(channelMatch{1}{1});
            col = str2double(channelMatch{1}{2});

            % Calculate new row index
            newRow = row + channelOffset;
            newName = sprintf('channels_%d_%d', newRow, col);

            % Add to result table
            result.(newName) = allData{i}.(oldName)(1);
        end

        % Update offset for next file
        if ~isempty(channelRowIndices)
            channelOffset = channelOffset + max(channelRowIndices);
        end
    end

    % 4. Process all coordinate columns next (in order)
    coordOffset = 0;
    for i = 1:length(allData)
        varNames = allData{i}.Properties.VariableNames;

        % Find all coord columns in this file
        coordCols = {};
        coordRowIndices = [];
        coordColIndices = [];

        for j = 1:length(varNames)
            varName = varNames{j};
            coordMatch = regexp(varName, coordColPattern, 'tokens');

            if ~isempty(coordMatch)
                row = str2double(coordMatch{1}{1});
                col = str2double(coordMatch{1}{2});

                coordCols{end+1} = varName;
                coordRowIndices(end+1) = row;
                coordColIndices(end+1) = col;
            end
        end

        % Sort the coord columns by row, then by column
        [~, sortIdx] = sortrows([coordRowIndices', coordColIndices']);
        coordCols = coordCols(sortIdx);

        % Add coord columns to result, renaming as needed
        for j = 1:length(coordCols)
            oldName = coordCols{j};
            coordMatch = regexp(oldName, coordColPattern, 'tokens');
            row = str2double(coordMatch{1}{1});
            col = str2double(coordMatch{1}{2});

            % Calculate new row index
            newRow = row + coordOffset;
            newName = sprintf('coords_%d_%d', newRow, col);

            % Add to result table
            result.(newName) = allData{i}.(oldName)(1);
        end

        % Update offset for next file
        if ~isempty(coordRowIndices)
            coordOffset = coordOffset + max(coordRowIndices);
        end
    end

    % 5. Finally, add all non-special columns
    for i = 1:length(allNonSpecialCols)
        colName = allNonSpecialCols{i};

        % Special handling for DivNm columns
        if strcmp(colName, 'DivNm')
            % First file's DivNm becomes DivNm_1 (with underscore)
            result.DivNm_1 = allData{1}.DivNm(1);

            % Add DivNm from other files
            for j = 2:length(allData)
                if ismember('DivNm', allData{j}.Properties.VariableNames)
                    newColName = sprintf('DivNm_%d', j);
                    result.(newColName) = allData{j}.DivNm(1);
                end
            end
        else
            % For all other columns, take value from first file that has it
            for j = 1:length(allData)
                if ismember(colName, allData{j}.Properties.VariableNames)
                    result.(colName) = allData{j}.(colName)(1);
                    break;
                end
            end
        end
    end

    % 6. Update specific fields again to ensure they're correct
    if ismember('outputDataFolderName', result.Properties.VariableNames)
        result.outputDataFolderName(1) = string(outputFolderName);
    end

    % Keep rawData from first CSV (don't modify it)
    % if ismember('rawData', result.Properties.VariableNames)
    %     result.rawData(1) = "default_value";
    % end

    if ismember('spreadSheetFileName', result.Properties.VariableNames)
        result.spreadSheetFileName(1) = string(divMergedCsvPath);
    end

    % Apply specialized handling to remove NaN values from all columns
    result = removeNaNValues(result);

    % Write the final merged data to output file with full precision
    % Using 'Precision', 'full' ensures all decimal places are preserved
    writetableOptions = {'Delimiter', ',', 'QuoteStrings', true, 'FileType', 'text'};
    
    % Use custom writer to preserve full decimal precision
    outputCsvFile = fullfile(outputFolder, ['Parameters_' outputFolderName '.csv']);
    writeTableWithFullPrecision(result, outputCsvFile);
    fprintf('Successfully wrote merged parameter CSV file to: %s\n', outputCsvFile);
end

function resultTable = removeNaNValues(inputTable)
    % removeNaNValues - Helper function to replace NaN values with appropriate empty values
    %
    % This function processes a table and replaces NaN values according to column type:
    % - For string columns: NaN becomes empty string ""
    % - For numeric columns: NaN is preserved but will be written as empty to CSV
    % - For cell columns: NaN becomes empty cell {}
    % - For categorical columns: NaN becomes <undefined>
    
    resultTable = inputTable;
    varNames = resultTable.Properties.VariableNames;
    
    % Create a temporary table to use during the transformation
    tempTable = table();
    
    % Process each column
    for i = 1:length(varNames)
        varName = varNames{i};
        columnData = resultTable.(varName);
        
        % Handle different data types
        if isa(columnData, 'string')
            % For string arrays, replace missing values with empty strings
            columnData(ismissing(columnData)) = "";
            tempTable.(varName) = columnData;
        elseif isa(columnData, 'cell')
            % For cell arrays, replace NaN cells with empty cells
            nanIdx = cellfun(@(x) isa(x, 'double') && isnan(x), columnData);
            if any(nanIdx)
                columnData(nanIdx) = {''};
            end
            tempTable.(varName) = columnData;
        elseif isa(columnData, 'double') || isa(columnData, 'single')
            % For numeric data, keep full precision but convert NaNs to empty
            if any(isnan(columnData))
                % Create a cell array to handle mixed numeric and empty data
                cellData = cell(size(columnData));
                for j = 1:numel(columnData)
                    if isnan(columnData(j))
                        cellData{j} = '';
                    else
                        cellData{j} = columnData(j);
                    end
                end
                tempTable.(varName) = cellData;
            else
                % No NaNs, keep as numeric for full precision
                tempTable.(varName) = columnData;
            end
        elseif isa(columnData, 'categorical')
            % For categorical, convert NaN to <undefined>
            columnData(isundefined(columnData)) = categorical({''});
            tempTable.(varName) = columnData;
        else
            % For other types, keep as is
            tempTable.(varName) = columnData;
        end
    end
    
    resultTable = tempTable;
end

function writeTableWithFullPrecision(data, filename)
    % Custom function to write table with full numeric precision
    
    % Get variable names
    varNames = data.Properties.VariableNames;
    numVars = length(varNames);
    
    % Open file for writing
    fileID = fopen(filename, 'w');
    if fileID == -1
        error('Could not open file for writing: %s', filename);
    end
    
    % Write header
    fprintf(fileID, '%s', varNames{1});
    for i = 2:numVars
        fprintf(fileID, ',%s', varNames{i});
    end
    fprintf(fileID, '\n');
    
    % Write data rows
    numRows = height(data);
    for row = 1:numRows
        % Process each column in the row
        for col = 1:numVars
            if col > 1
                fprintf(fileID, ',');
            end
            
            % Get the current value
            value = data{row, col};
            
            % Handle different types of data
            if iscell(value)
                % Cell value
                cellValue = value{1};
                if isempty(cellValue)
                    % Empty - write nothing
                    fprintf(fileID, '');
                elseif ischar(cellValue) || isstring(cellValue)
                    % Text - quote if needed
                    fprintf(fileID, '"%s"', strrep(char(cellValue), '"', '""'));
                elseif isnumeric(cellValue)
                    % Numeric - write with full precision
                    if isnan(cellValue)
                        fprintf(fileID, '');
                    else
                        fprintf(fileID, '%.15g', cellValue);
                    end
                else
                    % Other types
                    fprintf(fileID, '"%s"', char(string(cellValue)));
                end
            elseif ischar(value) || isstring(value)
                % Text - quote if needed
                if isempty(char(value))
                    fprintf(fileID, '');
                else
                    fprintf(fileID, '"%s"', strrep(char(value), '"', '""'));
                end
            elseif isnumeric(value)
                % Numeric - write with full precision
                if isnan(value)
                    fprintf(fileID, '');
                else
                    fprintf(fileID, '%.15g', value);
                end
            elseif iscategorical(value)
                % Categorical
                if isundefined(value)
                    fprintf(fileID, '');
                else
                    fprintf(fileID, '"%s"', char(value));
                end
            else
                % Other types
                fprintf(fileID, '"%s"', char(string(value)));
            end
        end
        fprintf(fileID, '\n');
    end
    
    % Close the file
    fclose(fileID);
end