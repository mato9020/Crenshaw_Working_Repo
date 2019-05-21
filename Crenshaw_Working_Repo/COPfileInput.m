classdef COPfileInput < matlab.System
    %COPfileInput is simply meant to split the .tsv or .forces file into
    %the raw data used for processing. This will automatically detect the
    %file input as long as it is .tsv or .forces. 

    properties
        filename;
        rawData;
        COPData;
        sRate;
    end
    
    methods
        function obj = COPfileInput(arg)
            obj.filename = (arg);
            [~,~,ext] = fileparts(obj.filename);
            if strcmp(ext,'.tsv') == 1
                obj.COPData = importdata(obj.filename,'\t',27);%Hardcoded for a specific TSV output (only force data and headers should be included)
            else
                obj.rawData = importdata(obj.filename,'\t');
                obj.COPData=obj.rawData.data(:,5:6);
            end
        end
        function setrawData(obj,value)
            obj.rawData = value;
        end
        function val = getCOPData(obj)
            val = obj.COPData;
        end
    end
end
