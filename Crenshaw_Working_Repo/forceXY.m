classdef forceXY < matlab.System
    %Object forceXY holds the data set and normalizes and filters the data
    %set.
    %The each data set from raw->unfilt->forceSet holds its value so it
    %does not get normalized or filtered more than once.
        
    properties(Hidden)
        COfreq = 10;%Cutpff frequency
        downsampleRate;
        downsampleFact = 2;
        sRate = 1200;%sampling Rate
    end
    properties
        
        endTime = 30;
        unfiltForceSet;%raw data normalized to mean value
        forceSet;%filtered forceplate data using butterworth filt
        %Prietto Vars%%%%%%%%%%%%%%%%%%%%%%
        RMS;
        MDIST;
        TOTEX;
        MVELO;
        MFREQ;
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    properties(Access = protected)
       rawForceSet;%raw data with coords relative to FP origin
    end
    properties(Logical,Access = private)%logical prop to identify when the force set is a transverse
        isHypot = false;
        isLowSrate = false;
    end
    methods
        function obj = forceXY(array,sRate,hypotArray)
            if nargin == 0
               obj.rawForceSet = [];
               obj.forceSet = [];
               obj.unfiltForceSet = [];
            elseif nargin ==2%%Filters and normalizes AP and ML arrays
                obj.rawForceSet = array;
                if(iscell(sRate))
                    sRate = str2num(cell2mat(sRate));
                end
                filtForceSet(obj,sRate);
                obj.sRate = sRate;
                set(obj,'unfiltForceSet',(obj.unfiltForceSet-mean(obj.unfiltForceSet)));
                set(obj,'forceSet',(obj.forceSet-mean(obj.forceSet)));
                calcVar(obj);%initializes all prietto vars for data array
            else
                obj.forceSet = hypot(array,hypotArray);%For construction of the transverse plane data array and calulations
                obj.isHypot = true;
                obj.TOTEX = sum( (diff(array).^2+diff(hypotArray).^2).^(1/2) );%calculates the totex of the transverse as per Prietto
                calcVar(obj);
            end
            
            
        end
        function calcVar(obj)%Base Calculations of the vars that only require the force data to calculate
            obj.MDIST = sum(abs(obj.forceSet))/length(obj.forceSet);
            obj.RMS = rms(obj.forceSet);
            if obj.isHypot == false%TOTEX calc for the TV is based on the ML and AP force data and thus is not calculated here
                obj.TOTEX = sum(abs(diff(obj.forceSet)));
            end
            obj.MVELO = obj.TOTEX/obj.endTime;
            obj.MFREQ = obj.MVELO/(2*pi*obj.MDIST);
        end
        function [RMS,MDIST,TOTEX,MVELO,MFREQ] = getForceVar(obj)%for base vars
            MDIST = get(obj,'MDIST');
            RMS = get(obj,'RMS');
            TOTEX = get(obj,'TOTEX');
            MVELO = get(obj,'MVELO');
            MFREQ = get(obj,'MFREQ');
            
        end
        function val = get.forceSet(obj)
            val = obj.forceSet;
        end
        function filtForceSet(obj,sRate)%function to filter the force set
            %downsample to 600hz
            if sRate == 1200
                
                aVal = floor(size(obj.rawForceSet,1)/(sRate*5));
                obj.unfiltForceSet = downsample(obj.rawForceSet(1:sRate*aVal*5,:),obj.downsampleFact);
            elseif sRate == 1000
                aVal = floor(size(obj.rawForceSet,1)/(sRate*5));
                obj.unfiltForceSet = resample(obj.rawForceSet(1:sRate*aVal*5,:),3,5);
            end
            %Filtering
            obj.endTime = aVal*5;
            obj.downsampleRate = length(obj.unfiltForceSet)/obj.endTime;%sRate/obj.downsampleFact;
            [z,p] = butter(2,obj.COfreq/(obj.downsampleRate/2));  %Assumes sampling frequency is 600 Hz, 2nd order, or forth order recursive
            obj.forceSet=filtfilt(z,p,obj.unfiltForceSet);%Filtering signal
        end
        function val = std(obj)%the standard deviation of the AP and ML is equal to the RMS value. If manually calculated std~=rms but only for AP and ML
           val = obj.RMS; 
        end
        function val = stdTV(obj)%Prietto standard deviation of RD time series
           val = (obj.RMS^2-obj.MDIST^2)^(1/2);
        end
        function robj = changeDur(obj,endTime)%Function to change the duration of the force set. Will change the unfilt and filt force set, but leave the raw force set property null. 
           robj = forceXY();
           set(robj,'endTime',endTime);
           set(robj,'unfiltForceSet',obj.unfiltForceSet(1:endTime*(obj.downsampleRate)));
           set(robj,'forceSet',obj.forceSet(1:endTime*obj.downsampleRate));
           calcVar(robj);
        end
    end
end
