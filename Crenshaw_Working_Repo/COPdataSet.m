classdef COPdataSet < matlab.System
    % Raw data path for FileInput object
    %This object holds the two x and y data sets and creates the transverse
    %data set "Hypot"
    %AP and ML are the X and Y coordinate sets respectively,Hypot is the
    %data set calculated as the length of the hypotenuse formed from the x
    %and y
%     

    % Public, tunable properties
    properties %
        AP;
        ML;
        Hypot;
        varTable=cell(3,12);
        endTime;
        subject;
    end


    methods
        function obj = COPdataSet(AParray,sRate,MLarray)%Constructor func for data set.
            if nargin == 0
               obj.AP= forceXY();
               obj.ML = forceXY();
               obj.Hypot = forceXY();
            else
                obj.AP = forceXY(AParray,sRate);
                obj.ML = forceXY(MLarray,sRate);
                obj.Hypot = forceXY(get(obj.AP,'forceSet'),sRate,get(obj.ML,'forceSet'));
            end
        end
        
        function table = blockBuild(obj,subject,condition,endTime)%Constructs the table used in the app figure
            AREACC=0;
            AREACE=0;
            AREASW=0;
            for i=1:3%loops through each object and calculates necessary vars
                [force,name] = dLoop(obj,i);
                [RMS,MDIST,TOTEX,MVELO,MFREQ] = getForceVar(force);
                if i == 3
                    AREASW = SwayArea(obj);
                    AREACC =ccArea(obj);
                    AREACE = ceArea(obj);
                    TOTEX = get(force,'TOTEX');
                end
                obj.varTable(i,:) = {subject,condition,endTime,name,RMS,MDIST,TOTEX,MVELO,MFREQ,AREACC,AREACE,AREASW};%Holds the value for each variable
                table = obj.varTable;
            end
        end

        function [returnObj,name] = dLoop(obj,i)%structure used to iterate through each object
            switch i
                case 1
                    name = 'AP';
                    returnObj = obj.AP;
                case 2 
                    name = 'ML';
                    returnObj = obj.ML;
                case 3
                    name = 'Hypot';
                    returnObj = obj.Hypot;
            end
        end
        %%%%%%%%%%%%%%%%%% Further Prietto Vars%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function area = ccArea(obj)
            rad = get(obj.Hypot,'MDIST')+1.645*stdTV(obj.Hypot);
            area = pi*(rad)^2;
        end
        function area = ceArea(obj)
            area = 6*pi*( (std(obj.AP)^2).*(std(obj.ML)^2)-((1/(length(get(obj.AP,'forceSet'))))*sum( get(obj.AP,'forceSet').*get(obj.ML,'forceSet') ))^2) ^(1/2);
        end
        function area = SwayArea(obj)
           ap = get(obj.AP,'forceSet');
           ml = get(obj.ML,'forceSet');
           accum = 0;
           for n =1:length(ap)-1
               accum = accum+abs((ap(n+1)*ml(n)-ml(n+1)*ap(n)));
           end
           area = (1/(2*get(obj.AP,'endTime')))*accum;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function robj = changeDur(obj,endTime)%general function to call the changeDur function of each forceXY in the data set. Creates a new data set and takes the amount of time requested
            robj = COPdataSet();
            for i=1:2%loops through each object and calculates necessary vars
                [~,name] = dLoop(obj,i);
                set(robj,name,changeDur(get(obj,name),endTime));
            end
            set(robj,'Hypot',(forceXY(get(robj.AP,'forceSet'),get(obj.ML,'sRate'),get(robj.ML,'forceSet'))));
        end
    end
end
