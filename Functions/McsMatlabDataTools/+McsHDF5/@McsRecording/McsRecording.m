classdef McsRecording < handle
% Stores a single recording.
%
% The different streams present in the recording are sorted into the
% {Analog,Frame,Event,Segment,TimeStamp}Stream fields where they are stored
% as cell arrays.
%
% (c) 2016 by Multi Channel Systems MCS GmbH
    
    properties (SetAccess = private)
        RecordingID = 0         % (scalar) The recording ID
        RecordingType           % (string) The recording type
        TimeStamp               % (scalar) Time stamp of the start date of the recording, in 10^{-7} s
        Duration                % (scalar) Duration of the recording, in 10^{-7} s
        Label                   % (string) The recording label
        Comment                 % (string) Comment
        AnalogStream = {};      % (cell array) McsAnalogStream objects, one for each analog stream
        FrameStream = {};       % (cell array) McsFrameStream objects, one for each frame stream
        EventStream = {};       % (cell array) McsEventStream objects, one for each event stream
        SegmentStream = {};     % (cell array) McsSegmentStream objects, one for each segment stream
        TimeStampStream = {};   % (cell array) McsTimeStampStream objects, one for each time stamp stream
    end
    
    methods
        
        function rec = McsRecording(filename, recStruct, cfg)
        % Reads a single recording inside a HDF5 file.
        %
        % function rec = McsRecording(filename, recStruct, cfg)
        %
        % Input:
        %   filename    -   (string) Name of the HDF5 file
        %   recStruct   -   The recording subtree of the structure 
        %                   generated by the h5info command
        %   cfg     -   (optional) configuration structure, contains one or
        %               more of the following fields:
        %               'dataType': The type of the data, can be one of
        %               'double' (default), 'single' or 'raw'. For 'double'
        %               and 'single' the data is converted to meaningful
        %               units, while for 'raw' no conversion is done and
        %               the data is kept in ADC units. This uses less
        %               memory than the conversion to double, but you might
        %               have to convert the data prior to analysis, for
        %               example by using the getConvertedData function.
        %               'timeStampDataType': The type of the time stamps,
        %               can be either 'int64' (default) or 'double'. Using
        %               'double' is useful for older Matlab version without
        %               int64 arithmetic.
        %               'correctConversionFactorOrientation': (bool) flag
        %               that signals whether to transpose the orientation
        %               of the conversion factor matrix in frame streams.
        %               This is necessary for older DataManager versions
        %
        % Output:
        %   rec         -   A McsRecording object
        %
            if exist('h5info')
                mode = 'h5';
            else
                mode = 'hdf5';
            end
            source = 'DataManager';
            dataAttributes = recStruct.Attributes;
            for fni = 1:length(dataAttributes)
                if strcmp(mode,'h5')
                    rec.(dataAttributes(fni).Name) = dataAttributes(fni).Value;
                else
                    str = regexp(dataAttributes(fni).Name,'/\w+$','match');
                    if isa(dataAttributes(fni).Value,'hdf5.h5string')
                        rec.(str{length(str)}(2:end)) = dataAttributes(fni).Value.Data;
                    else
                        rec.(str{length(str)}(2:end)) = dataAttributes(fni).Value;
                    end
                end
            end

            for gidx = 1:length(recStruct.Groups)
                groupname = recStruct.Groups(gidx).Name;
                
                if ~isempty(strfind(groupname,'AnalogStream'))
                    % read analog streams
                    count = 1;
                    for streams = 1:length(recStruct.Groups(gidx).Groups)
                        if length(recStruct.Groups(gidx).Groups(streams).Datasets) <= 1
                            continue;
                        end
                        rec.AnalogStream{count} = McsHDF5.McsAnalogStream(filename, recStruct.Groups(gidx).Groups(streams), source, cfg);
                        count = count + 1;
                    end
                    
                elseif ~isempty(strfind(groupname,'FrameStream'))
                    % read frame streams
                    % can't do the data set count here because the frame
                    % data entities are below the frame stream
                    for streams = 1:length(recStruct.Groups(gidx).Groups)
                        rec.FrameStream{streams} = McsHDF5.McsFrameStream(filename, recStruct.Groups(gidx).Groups(streams), source, cfg);
                    end
                    
                elseif ~isempty(strfind(groupname,'EventStream'))
                    % read event streams
                    count = 1;
                    for streams = 1:length(recStruct.Groups(gidx).Groups)
                        if length(recStruct.Groups(gidx).Groups(streams).Datasets) <= 1
                            continue;
                        end
                        rec.EventStream{count} = McsHDF5.McsEventStream(filename, recStruct.Groups(gidx).Groups(streams), source, cfg);
                        count = count + 1;
                    end
                    
                elseif ~isempty(strfind(groupname,'SegmentStream'))
                    % read segment streams
                    count = 1;
                    for streams = 1:length(recStruct.Groups(gidx).Groups)
                        if length(recStruct.Groups(gidx).Groups(streams).Datasets) <= 1
                            continue;
                        end
                        rec.SegmentStream{count} = McsHDF5.McsSegmentStream.makeSegmentStream(filename, recStruct.Groups(gidx).Groups(streams), cfg);
                        count = count + 1;
                    end
                
                elseif ~isempty(strfind(groupname,'TimeStampStream'))
                    % read timestamp streams
                    count = 1;
                    for streams = 1:length(recStruct.Groups(gidx).Groups)
                        if length(recStruct.Groups(gidx).Groups(streams).Datasets) <= 1
                            continue;
                        end
                        rec.TimeStampStream{count} = McsHDF5.McsTimeStampStream(filename, recStruct.Groups(gidx).Groups(streams), cfg);
                        count = count + 1;
                    end
                end      
            end    
        end
    end
end