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