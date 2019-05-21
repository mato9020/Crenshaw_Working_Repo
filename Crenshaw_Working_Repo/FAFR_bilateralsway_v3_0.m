%Program to analyze and output COP variables of postural sway for the FAFR
%study
%v1-0 Jeremy Crenshaw 4/29/2013
%v2-0 Jeremy Crenshaw 5/17/2013 to use updated COP and COM calculations from
%newer Vis3d pipelines.  COM calculated from markers to include head and upper extremities
%v3-0 Jeremy Crenshaw 7/2/2013 to implement stochastic variables as calculated by programs by V Lugade



clear all;
close all;

%TO DO:
%Keep an eye on total power and U0... not close to prieto 1996 values.

%Change the first subject and last subject of range to analyze
firstsubject=124;
lastsubject=125;

newoutput=[];

%Output data options y=yes n=no (other inputs will be no)
fileoutput='y';

%Cut-off freq for lowpass filter of COP and GRF
COfreq=5;
time=30;

%Convert all to mm? 'y' for yes (no is  in m)
convert='y';

if convert=='y';
    cf=1000; %conversion factor
else
    cf=1;
end

%All text files are in
%I:\INVIVO\RESEARCH\FAFR_11-006984_Amin\MOTIONVIDEO\Processed-Data\fafra###a\
%Need following text files exported using Vis3D:
%fafra###a_sway_uni_MARKERS.txt
%fafra###a_sway_uni_OUTPUT.txt
%fafra###a_sway_open1_MARKERS.txt
%fafra###a_sway_open1_OUTPUT.txt
%fafra###a_sway_closed1_MARKERS.txt
%fafra###a_sway_closed1_OUTPUT.txt
%patientinfo.xls (NOT NEEDED FOR NOW)
%userInputs.csv
%Assumes all data is sampled at 120 Hz.
%needs calc_vel_acc.m, sgsdf.m routines in same folder; not for now, but if COMpos and vel are too noisy
%need MyPSD.m
%need com_calc_FAFR.m   same as com_calc.m by VA Lugade
%X is mediolateral, Y is anteroposterior.  subject faces +Y, -X is to subject's left.

%folder with all data
bigpath='I:\INVIVO\RESEARCH\FAFR_11-006984_Amin\MOTIONVIDEO\Processed-Data_Cortex4ONLY\';
outfile=['I:\INVIVO\RESEARCH\FAFR_11-006984_Amin\MOTIONVIDEO\CompiledData\Sway\','FAFRa_SwayOutput.xls'];
outfolder='I:\INVIVO\RESEARCH\FAFR_11-006984_Amin\MOTIONVIDEO\CompiledData\Sway\';




%Filtering considerations:
%Berg/Maki 1992, COP low-pass filtered (-3dB at 10HZ) and sampled at a rate of 50 Hz, same as Maki 1991)
%Chiari 2002, moments and forces filtered at 8 hz (30th order low-pass FIR filter with zero-phase) and down-sampled at 20 Hz
%Lafond 2004, force plafgorm signals filtered with a zero-lag 6th order butterworth filter at 10hz.
%Lin 2008, forces and moments low-pass filtered (2nd order, zero-phase-lag, butterworth, 5 Hz cut-off frequency)
%Maurer 2005, COP trace filtered with a first-order low=pass filter with a cut-off frequency of 5 Hz.

%continue with 5 Hz cut-off of COP data.  Record freq content under 5 Hz.


%loop to go through each subject's folder
 for i = firstsubject:lastsubject
%for i = [67]
    
    clearvars -except i firstsubject lastsubject fileoutput COfreq time convert cf bigpath newoutput outfile outfolder;
    close all;
    %Identify files needed
    filepath=[bigpath,'fafra',num2str(i,'%03d'),'a\'];
    CLOSEDmarkerfile=[filepath,'fafra',num2str(i,'%03d'),'a_sway_closed1_MARKERS.txt'];
    CLOSEDoutputfile=[filepath,'fafra',num2str(i,'%03d'),'a_sway_closed1_OUTPUT.txt'];
    OPENmarkerfile=[filepath,'fafra',num2str(i,'%03d'),'a_sway_open1_MARKERS.txt'];
    OPENoutputfile=[filepath,'fafra',num2str(i,'%03d'),'a_sway_open1_OUTPUT.txt'];
    UNImarkerfile=[filepath,'fafra',num2str(i,'%03d'),'a_sway_uni_MARKERS.txt'];
    UNIoutputfile=[filepath,'fafra',num2str(i,'%03d'),'a_sway_uni_OUTPUT.txt'];
    patientinfofile=[filepath,'patientinfo.xls'];
    userinputfile=[filepath,'userInputs.csv'];
    
    %NOTE: with switch to cortex, closed01 and open01... make program flexible.
    %check to see if files exist, if they don't change number from "1" ot "01"
    if exist(CLOSEDmarkerfile,'file')==0 %does not exist, change to 01
        CLOSEDmarkerfile=[filepath,'fafra',num2str(i,'%03d'),'a_sway_closed01_MARKERS.txt'];
        CLOSEDoutputfile=[filepath,'fafra',num2str(i,'%03d'),'a_sway_closed01_OUTPUT.txt'];
        OPENmarkerfile=[filepath,'fafra',num2str(i,'%03d'),'a_sway_open01_MARKERS.txt'];
        OPENoutputfile=[filepath,'fafra',num2str(i,'%03d'),'a_sway_open01_OUTPUT.txt'];
    end
    %disp(CLOSEDmarkerfile); %troublshooting line
    
%     if i==78
%         CLOSEDmarkerfile=[filepath,'fafra',num2str(i,'%03d'),'a_swayClosed02_MARKERS.txt'];
%         CLOSEDoutputfile=[filepath,'fafra',num2str(i,'%03d'),'a_swayClosed02_OUTPUT.txt'];
%     end
    
    if i==86
        OPENmarkerfile=[filepath,'fafra',num2str(i,'%03d'),'a_sway_open02_MARKERS.txt'];
        OPENoutputfile=[filepath,'fafra',num2str(i,'%03d'),'a_sway_open02_OUTPUT.txt'];
    end
    
    subject=i;
    disp(subject)
    %Access user input file for foot width and length
    [UserInputNum,UserInputTxt,UserInputRaw]=xlsread(userinputfile);
    RFootL=UserInputNum(6)/100*cf; %m
    LFootL=UserInputNum(7)/100*cf; %m
    RFootW=UserInputNum(8)/100*cf; %m
    LFootW=UserInputNum(9)/100*cf; %m
    height=UserInputNum(2)/100*cf; %m
    Pheight=height*0.586;  %pendulum height from from Clauser 1969
    AvgfootL=mean([RFootL,LFootL]);
    AvgfootW=mean([RFootW,LFootW]);
    
    
    %Recreate BOS as per Chiari 2002 for bilateral sway
    %access open, closed sway marker files.
    OPENmarkers=importdata(OPENmarkerfile);
    OPENmarkersize=size(OPENmarkers.data);
    OPENmarkerheaders=textscan(OPENmarkers.textdata{2,1}, '%s', OPENmarkersize(2), 'delimiter', '\t');
    CLOSEDmarkers=importdata(CLOSEDmarkerfile);
    CLOSEDmarkersize=size(CLOSEDmarkers.data);
    CLOSEDmarkerheaders=textscan(CLOSEDmarkers.textdata{2,1}, '%s', CLOSEDmarkersize(2), 'delimiter', '\t');
    OPENLheel=OPENmarkers.data(:,14:16)*cf;
    OPENLtoe=OPENmarkers.data(:,56:58)*cf;
    OPENRheel=OPENmarkers.data(:,71:73)*cf;
    OPENRtoe=OPENmarkers.data(:,116:118)*cf;
    CLOSEDLheel=CLOSEDmarkers.data(:,14:16)*cf;
    CLOSEDLtoe=CLOSEDmarkers.data(:,56:58)*cf;
    CLOSEDRheel=CLOSEDmarkers.data(:,71:73)*cf;
    CLOSEDRtoe=CLOSEDmarkers.data(:,116:118)*cf;
    
    %Create vector from heel to toe
    OPENRfootv=mean(OPENRtoe-OPENRheel);
    OPENLfootv=mean(OPENLtoe-OPENLheel);
    CLOSEDRfootv=mean(CLOSEDRtoe-CLOSEDRheel);
    CLOSEDLfootv=mean(CLOSEDLtoe-CLOSEDLheel);
    
       
    %mean heel position
    OPENRheelmean=mean(OPENRheel(:,1:2));
    OPENLheelmean=mean(OPENLheel(:,1:2));
    CLOSEDRheelmean=mean(CLOSEDRheel(:,1:2));
    CLOSEDLheelmean=mean(CLOSEDLheel(:,1:2));
    
    %Magnitude of this vector
    MagOPENRfootv=sqrt(OPENRfootv(1).^2+OPENRfootv(2).^2);
    MagOPENLfootv=sqrt(OPENLfootv(1).^2+OPENLfootv(2).^2);
    MagCLOSEDRfootv=sqrt(CLOSEDRfootv(1).^2+CLOSEDRfootv(2).^2);
    MagCLOSEDLfootv=sqrt(CLOSEDLfootv(1).^2+CLOSEDLfootv(2).^2);
    
    %foot unit vectors
    OPENRfootu=OPENRfootv(1:2)./MagOPENRfootv;
    OPENLfootu=OPENLfootv(1:2)./MagOPENLfootv;
    CLOSEDRfootu=CLOSEDRfootv(1:2)./MagCLOSEDRfootv;
    CLOSEDLfootu=CLOSEDLfootv(1:2)./MagCLOSEDLfootv;
    
    
    
    %Create virtual "big toe" markers (will be more in middle of foot)
    OPENRBT=OPENRheelmean+RFootL*OPENRfootu;
    OPENLBT=OPENLheelmean+LFootL*OPENLfootu;
    CLOSEDRBT=CLOSEDRheelmean+RFootL*CLOSEDRfootu;
    CLOSEDLBT=CLOSEDLheelmean+LFootL*CLOSEDLfootu;
    
    OPENbos=[OPENRBT; OPENLBT; OPENLheelmean; OPENRheelmean; OPENRBT];
    CLOSEDbos=[CLOSEDRBT; CLOSEDLBT; CLOSEDLheelmean; CLOSEDRheelmean; CLOSEDRBT];
    
    %Vectors for beta (angle between vector from left to right heel and x axis)
    OPENheelV=(OPENRheelmean-OPENLheelmean);
    OPENheelmag=sqrt(OPENheelV(1).^2+OPENheelV(2).^2);
    OPENheelU=OPENheelV./OPENheelmag;
    CLOSEDheelV=(CLOSEDRheelmean-CLOSEDLheelmean);
    CLOSEDheelmag=sqrt(CLOSEDheelV(1).^2+CLOSEDheelV(2).^2);
    CLOSEDheelU=CLOSEDheelV./CLOSEDheelmag;
        
    %As modified from Chieari 2002, calculate alpha (foot duck angle)
    % calculate BOS length, BOS width, and BOS area
    OPENfootalpha=rad2deg(acos(dot(OPENRfootu,OPENLfootu)));
    CLOSEDfootalpha=rad2deg(acos(dot(CLOSEDRfootu,CLOSEDLfootu)));
    OPENfootbeta=rad2deg(acos(dot(OPENheelU,[1 0])));
    CLOSEDfootbeta=rad2deg(acos(dot(CLOSEDheelU,[1 0])));    
    OPENbosL=max(OPENbos(:,2))-min(OPENbos(:,2));
    OPENbosW=max(OPENbos(:,1))-min(OPENbos(:,1));
    CLOSEDbosL=max(CLOSEDbos(:,2))-min(CLOSEDbos(:,2));
    CLOSEDbosW=max(CLOSEDbos(:,1))-min(CLOSEDbos(:,1));
    %Area of BOS using polyarea
    OPENbosA=polyarea(OPENbos(1:4,1),OPENbos(1:4,2));
    CLOSEDbosA=polyarea(CLOSEDbos(1:4,1),CLOSEDbos(1:4,2));
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Access uni, open, closed sway output files.  UNIoutput.data will have
    %numbers
    UNIoutput=importdata(UNIoutputfile);
    OPENoutput=importdata(OPENoutputfile);
    CLOSEDoutput=importdata(CLOSEDoutputfile);
    
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get Kicking Limb out of Unilateral output file
    KickLimb=UNIoutput.data(1,2); %1 is Right, -1 is Left
    
    %Get max L and R side unilateral stance time, make relative to kicking,
    %non-kicking limb
    %Determine how many trials were taken
    UNIoutputsize=size(UNIoutput.data);
    
    UNItimes=zeros(floor(UNIoutputsize(2)/15),2);
    for j = 1:floor(UNIoutputsize(2)/15)  %should elicit 1:# of trials
        UNItimes(j,1)=UNIoutput.data(1,(3+15*(j-1))); %stance limb 1=R -1=L
        %gets last time of time column, data collected until foot goes down or subject "hops", then subtracts beginning time
        UNItimes(j,2)=UNIoutput.data(sum(isfinite(UNIoutput.data(:,4+15*(j-1)))),4+15*(j-1))-UNIoutput.data(1,4+15*(j-1));
        %times in s
    end
    %Determine max time for all, max time for kicking limb, max time for non-kicking limb, and first stance limb
    FirstUNIstance=UNItimes(1,1);
    MaxUNIboth=max(UNItimes(:,2));
    MaxUNIright=max((UNItimes(:,1)==1).*(UNItimes(:,2)));
    MaxUNIleft=max((UNItimes(:,1)==-1).*(UNItimes(:,2)));
    if KickLimb==1 %Right is the kicking limb
        MaxUNIkick=MaxUNIright;
        MaxUNInonkick=MaxUNIleft;
    else %Left is the kicking limb
        MaxUNIkick=MaxUNIleft;
        MaxUNInonkick=MaxUNIright;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %mean GRF asymmetry (kicking/non-kicking)
    OPENGRFasymm=mean(OPENoutput.data(:,11));
    CLOSEDGRFasymm=mean(CLOSEDoutput.data(:,11));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %get collection frequency (should be 120, but safegaurd)
    OPENfreq=round(1/(OPENoutput.data(2,2)-OPENoutput.data(1,2)));
    CLOSEDfreq=round(1/(CLOSEDoutput.data(2,2)-CLOSEDoutput.data(1,2)));
    
    if OPENfreq~=120 || CLOSEDfreq~=120
        disp('Frequency not 120 for a collection!')
        disp('Subject Number')
        disp(subject)
        disp('Eyes open collection frequency')
        disp(OPENfreq)
        disp('Eyes closed collection frequency')
        disp(CLOSEDfreq)
        disp('Proceeding with displayed collection frequencies');
    end
        
    
    %need to get COPnet global location to calculate pendulum length for xCOM calculation
    OPENCOPnetposG=OPENoutput.data(:,7:8)*cf;
    CLOSEDCOPnetposG=CLOSEDoutput.data(:,7:8)*cf;
    
    %for some trials, these variables have NaN...doesn't work if no NaN
%     EndFrameOPEN=find(isnan(OPENCOPnetposG),1,'first')-1;
%     OPENCOPnetposG=OPENCOPnetposG(1:EndFrameOPEN,:);
%     EndFrameCLOSED=find(isnan(CLOSEDCOPnetposG),1,'first')-1;
%     CLOSEDCOPnetposG=CLOSEDCOPnetposG(1:EndFrameCLOSED,:);
    
    %Trunk angles x(fore-aft pitch) y (left right roll), Gill 2001
    OPENtrunkAngle=OPENoutput.data(:,28:30);
    CLOSEDtrunkAngle=CLOSEDoutput.data(:,28:30);
    OPENtrunksway=max(OPENtrunkAngle)-min(OPENtrunkAngle);
    CLOSEDtrunksway=max(CLOSEDtrunkAngle)-min(CLOSEDtrunkAngle);
    OPENtrunkAngVel=diff(OPENtrunkAngle).*OPENfreq;
    CLOSEDtrunkAngVel=diff(CLOSEDtrunkAngle).*CLOSEDfreq;
    OPENtrunkswayvel=max(OPENtrunkAngVel)-min(OPENtrunkAngVel);
    CLOSEDtrunkswayvel=max(CLOSEDtrunkAngVel)-min(CLOSEDtrunkAngVel);
  
    

    
    %Frequency analysis of COP    
    %Low pass filter at 5 Hz, Maurer 2005... report %AUC at 5 Hz for, 1st order... but becomes 2nd order with filtfilt
    [Oz,Op] = butter(1,COfreq/(OPENfreq/2));
    [Cz,Cp] = butter(1,COfreq/(CLOSEDfreq/2)); %in case closed is collected at a different frequency
    for k=1:2 %x and y directions for COP
        [OPENPSD(:,k),OPENPSDfreq(:,k)] =MyPSD(OPENCOPnetposG(:,k),OPENfreq);
        OPENPSDauc=trapz(OPENPSDfreq(1:(COfreq+1),k),OPENPSD(1:(COfreq+1),k));
        OPENPSDtotalauc=trapz(OPENPSDfreq(:,k),OPENPSD(:,k));
        OPENPSDpercentauc(k)=(OPENPSDauc/OPENPSDtotalauc)*100;
        [CLOSEDPSD(:,k),CLOSEDPSDfreq(:,k)] =MyPSD(CLOSEDCOPnetposG(:,k),CLOSEDfreq);
        CLOSEDPSDauc=trapz(CLOSEDPSDfreq(1:(COfreq+1),k),CLOSEDPSD(1:(COfreq+1),k));
        CLOSEDPSDtotalauc=trapz(CLOSEDPSDfreq(:,k),CLOSEDPSD(:,k));
        CLOSEDPSDpercentauc(k)=(CLOSEDPSDauc/CLOSEDPSDtotalauc)*100;
        
        %Low-pass filter 5 hz, filtfilt doubles order from filter design of butter command above

        OPENCOPnetposGfilt(:,k)=filtfilt(Oz,Op,OPENCOPnetposG(:,k));
        CLOSEDCOPnetposGfilt(:,k)=filtfilt(Cz,Cp,CLOSEDCOPnetposG(:,k));
        
    end
    [OPENCOPpos_alt OPENCOPvel OPENCOPaccel]=calc_vel_acc(OPENCOPnetposGfilt,OPENfreq);
    [CLOSEDCOPpos_alt CLOSEDCOPvel CLOSEDCOPaccel]=calc_vel_acc(CLOSEDCOPnetposGfilt,CLOSEDfreq);
    
    %COM variables... modified to use all markers
%     OPENCOMpos=OPENoutput.data(:,13:15)*cf;
%     CLOSEDCOMpos=CLOSEDoutput.data(:,13:15)*cf;
    [OPENCOMpos OPENsegCOMpos] = com_calc_FAFR(OPENmarkers.data,OPENmarkerheaders{1});
    [CLOSEDCOMpos CLOSEDsegCOMpos] = com_calc_FAFR(CLOSEDmarkers.data,CLOSEDmarkerheaders{1});
    OPENCOMpos=OPENCOMpos*cf;
    OPENsegCOMpos=OPENsegCOMpos*cf;
    CLOSEDCOMpos=CLOSEDCOMpos*cf;
    CLOSEDsegCOMpos=CLOSEDsegCOMpos*cf;
    
    OPENCOMposUF=OPENCOMpos; %for comparison later if needed
    CLOSEDCOMposUF=CLOSEDCOMpos;
    for k=1:3 %x and y directions..and z... move this to before pend length if going with it
        [OPENCOMposPSD(:,k),OPENCOMposPSDfreq(:,k)] =MyPSD(OPENCOMpos(:,k),OPENfreq);
        OPENCOMposPSDauc=trapz(OPENCOMposPSDfreq(1:(COfreq+1),k),OPENCOMposPSD(1:(COfreq+1),k));
        OPENCOMposPSDtotalauc=trapz(OPENCOMposPSDfreq(:,k),OPENCOMposPSD(:,k));
        OPENCOMposPSDpercentauc(k)=(OPENCOMposPSDauc/OPENCOMposPSDtotalauc)*100;
        [CLOSEDCOMposPSD(:,k),CLOSEDCOMposPSDfreq(:,k)] =MyPSD(CLOSEDCOMpos(:,k),CLOSEDfreq);
        CLOSEDCOMposPSDauc=trapz(CLOSEDCOMposPSDfreq(1:(COfreq+1),k),CLOSEDCOMposPSD(1:(COfreq+1),k));
        CLOSEDCOMposPSDtotalauc=trapz(CLOSEDCOMposPSDfreq(:,k),CLOSEDCOMposPSD(:,k));
        CLOSEDCOMposPSDpercentauc(k)=(CLOSEDCOMposPSDauc/CLOSEDCOMposPSDtotalauc)*100;
        
        %Low-pass filter 5 hz, filtfilt doubles order from filter design of butter command above

        OPENCOMpos(:,k)=filtfilt(Oz,Op,OPENCOMpos(:,k));
        CLOSEDCOMpos(:,k)=filtfilt(Cz,Cp,CLOSEDCOMpos(:,k));
        
    end
    
    
    %calculate pendulum length, or distance from COP location to COM location

    OPENpendlength=sqrt((OPENCOMpos(:,1)-OPENCOPnetposG(:,1)).^2+(OPENCOMpos(:,2)-OPENCOPnetposG(:,2)).^2+(OPENCOMpos(:,3)).^2);
    CLOSEDpendlength=sqrt((CLOSEDCOMpos(:,1)-CLOSEDCOPnetposG(:,1)).^2+(CLOSEDCOMpos(:,2)-CLOSEDCOPnetposG(:,2)).^2+(CLOSEDCOMpos(:,3)).^2);
    meanOPENpendlength=mean(OPENpendlength);
    meanCLOSEDpendlength=mean(CLOSEDpendlength);
    
    %Calculate COM velocity using Savitzky-Golay algorithm (Savitzky and Golay, 1964) to low-pass filter velocity
%     OPENCOMvel=diff(OPENCOMpos).*OPENfreq; %for no filtering... but vel is too messy
%     CLOSEDCOMvel=diff(CLOSEDCOMpos).*CLOSEDfreq;
  
    [OPENCOMpos OPENCOMvel OPENCOMaccel]=calc_vel_acc(OPENCOMpos,OPENfreq);
    [CLOSEDCOMpos CLOSEDCOMvel CLOSEDCOMaccel]=calc_vel_acc(CLOSEDCOMpos,CLOSEDfreq);
    
    
 
   
    %calc_vel_acc
    %calculate extrapolated center of mass (xCOM) position (Hof, Condition of Dynamic Stability 2005)
    OPENxCOMpos(:,1)=OPENCOMpos(:,1)+OPENCOMvel(:,1).*sqrt(OPENpendlength/(9.80665*cf));
    OPENxCOMpos(:,2)=OPENCOMpos(:,2)+OPENCOMvel(:,2).*sqrt(OPENpendlength/(9.80665*cf));
    CLOSEDxCOMpos(:,1)=CLOSEDCOMpos(:,1)+CLOSEDCOMvel(:,1).*sqrt(CLOSEDpendlength/(9.80665*cf));
    CLOSEDxCOMpos(:,2)=CLOSEDCOMpos(:,2)+CLOSEDCOMvel(:,2).*sqrt(CLOSEDpendlength/(9.80665*cf));
    
    %calculate extrapolated center of pressure(xCOP) positition
    OPENxCOPpos(:,1)=OPENCOPnetposGfilt(:,1)+OPENCOPvel(:,1).*sqrt(Pheight/(9.80665*cf));
    OPENxCOPpos(:,2)=OPENCOPnetposGfilt(:,2)+OPENCOPvel(:,2).*sqrt(Pheight/(9.80665*cf));
    CLOSEDxCOPpos(:,1)=CLOSEDCOPnetposGfilt(:,1)+CLOSEDCOPvel(:,1).*sqrt(Pheight/(9.80665*cf));
    CLOSEDxCOPpos(:,2)=CLOSEDCOPnetposGfilt(:,2)+CLOSEDCOPvel(:,2).*sqrt(Pheight/(9.80665*cf));
    
    figure;
    plot(OPENxCOPpos(:,1),OPENxCOPpos(:,2),'b');
    title(['Sub ',num2str(subject),': Eyes Open: COP(black), xCOP (blue), COM (red), xCOM (green)']);
    hold on
    plot(OPENCOMpos(:,1),OPENCOMpos(:,2),'r')
    plot(OPENxCOMpos(:,1),OPENxCOMpos(:,2),'g')
    plot(OPENCOPnetposGfilt(:,1),OPENCOPnetposGfilt(:,2),'k')
    plot(OPENbos(:,1),OPENbos(:,2),'k--')
    
    figure;
    plot(CLOSEDxCOMpos(:,1),CLOSEDxCOMpos(:,2),'g')
    title(['Sub ',num2str(subject),': Eyes Closed: COP(black), xCOP (blue), COM (red), xCOM (green)']);
    hold on
    plot(CLOSEDCOMpos(:,1),CLOSEDCOMpos(:,2),'r')
    plot(CLOSEDxCOPpos(:,1),CLOSEDxCOPpos(:,2),'b');
    plot(CLOSEDCOPnetposGfilt(:,1),CLOSEDCOPnetposGfilt(:,2),'k')
    plot(CLOSEDbos(:,1),CLOSEDbos(:,2),'k--')

    
    %find mean of COP, COM, and xCOM, zero all time signals.
    for k=1:2
        OPENCOPzero(:,k)=OPENCOPnetposGfilt(:,k)-mean(OPENCOPnetposGfilt(:,k));
        CLOSEDCOPzero(:,k)=CLOSEDCOPnetposGfilt(:,k)-mean(CLOSEDCOPnetposGfilt(:,k));
        OPENCOMzero(:,k)=OPENCOMpos(:,k)-mean(OPENCOMpos(:,k));
        CLOSEDCOMzero(:,k)=CLOSEDCOMpos(:,k)-mean(CLOSEDCOMpos(:,k));
        OPENxCOMzero(:,k)=OPENxCOMpos(:,k)-mean(OPENxCOMpos(:,k));
        CLOSEDxCOMzero(:,k)=CLOSEDxCOMpos(:,k)-mean(CLOSEDxCOMpos(:,k));
        OPENxCOPzero(:,k)=OPENxCOPpos(:,k)-mean(OPENxCOPpos(:,k));
        CLOSEDxCOPzero(:,k)=CLOSEDxCOPpos(:,k)-mean(CLOSEDxCOPpos(:,k));
    end
    
    figure;
    plot(OPENxCOPzero(:,1),OPENxCOPzero(:,2),'b');
    title(['Sub ',num2str(subject),': Eyes Open: COP(black), xCOP (blue), COM (red), xCOM (green)']);
    hold on
    plot(OPENCOMzero(:,1),OPENCOMzero(:,2),'r')
    plot(OPENxCOMzero(:,1),OPENxCOMzero(:,2),'g')
    plot(OPENCOPzero(:,1),OPENCOPzero(:,2),'k')
    
    figure;
    plot(CLOSEDxCOMzero(:,1),CLOSEDxCOMzero(:,2),'g')
    title(['Sub ',num2str(subject),': Eyes Closed: COP(black), xCOP (blue), COM (red), xCOM (green)']);
    hold on
    plot(CLOSEDCOMzero(:,1),CLOSEDCOMzero(:,2),'r')
    plot(CLOSEDxCOPzero(:,1),CLOSEDxCOPzero(:,2),'b');
    plot(CLOSEDCOPzero(:,1),CLOSEDCOPzero(:,2),'k')
    
    freq=[OPENfreq OPENfreq OPENfreq OPENfreq CLOSEDfreq CLOSEDfreq CLOSEDfreq CLOSEDfreq];
    %Create matrix of all variables and conditions for easy loop iterations
    %Page Order, Eyes Open COP xCOP COM xCOM, then eyes closed COP xCOP COM xCOM, columns are x y
    paths=cat(3, OPENCOPzero, OPENxCOPzero, OPENCOMzero, OPENxCOMzero, CLOSEDCOPzero, CLOSEDxCOPzero, CLOSEDCOMzero, CLOSEDxCOMzero);
    vels=cat(3, OPENCOPvel, OPENCOMvel(:,1:2), CLOSEDCOPvel, CLOSEDCOMvel(:,1:2));
    for j=1:8 %type of signal
               
        for k=1:2
            %sway path length
            swaylength(j,k)=trapz(abs(diff(paths(:,k,j)))); %no x component needed for trapz because diff and trapz all rel to frames
            %mean difference (mean of absolute value)
            meandiff(j,k)=mean(abs(paths(:,k,j)));
            %root mean square distance (RMS)
            rmserror(j,k)=rms(paths(:,k,j));
            %max displacement
            maxdisp(j,k)=max(paths(:,k,j))-min(paths(:,k,j));
            %mean frequency (as per prieto_1996)
            meanfreq(j,k)=swaylength(j,k)/(4*sqrt(2)*time*meandiff(j,k));
            %power spectral density computed using the sinusoidal multitaper method with 8 tapers (as per prieto_1996)
            %Pdens is in units of power per radians per sample, Pfreq is in rad/sample.
            nw=(8+1)/2; %8 tapers
            [Pdens(:,k,j),Pfreq(:,k,j)]=pmtm(paths(:,k,j),nw,[],freq(j));
            %freq domain measures from 0.1465 to 5.01 Hz (Prieto does 0.15 to 5 Hz)
            %find indices closes to 0.15 hx and 5 Hx
            [MINlow MINlowi]=min(abs(Pfreq(:,k,j)-0.15));
            [MINhigh MINhighi]=min(abs(Pfreq(:,k,j)-5));
            TotPower(j,k)=trapz(Pfreq(MINlowi:MINhighi,k,j),Pdens(MINlowi:MINhighi,k,j));
            %keep in mind that units are in m, so values smaller than if were mm (multiply by 1000000?)
            %also, sampling frequency can affect value? higher sampling frequency produces higher values.
            
            %Find 50 and 95% power frequency
            [min50 min50i]=min(abs(cumtrapz(Pfreq(MINlowi:MINhighi,k,j),Pdens(MINlowi:MINhighi,k,j))./TotPower(j,k)-0.5));
            Power50(j,k)=Pfreq(MINlowi+min50i);
            [min95 min95i]=min(abs(cumtrapz(Pfreq(MINlowi:MINhighi,k,j),Pdens(MINlowi:MINhighi,k,j))./TotPower(j,k)-0.95));
            Power95(j,k)=Pfreq(MINlowi+min95i);
            
            %centroidal frequency and frequency dispersion as per Prieto, 1996 (equations 29 and 30)
            changeF=Pfreq(2)-Pfreq(1);
            hh=1;
            for h=MINlowi:MINhighi
                U2parts(hh)=(h*changeF)^2*Pdens(h,k,j);
                U1parts(hh)=(h*changeF)^1*Pdens(h,k,j);
                U0parts(hh)=Pdens(h,k,j);
                U0fparts(hh)=Pdens(h,k,j)*Pfreq(h,k,j);
                hh=hh+1;
            end
            U2=sum(U2parts);
            U1=sum(U1parts);
            U0=sum(U0parts);
            CFREQ(j,k)=sqrt(U2/U0);
            FREQD(j,k)=sqrt(1-(U1^2)/(U2*U0));
            U0total(j,k)=U0; 
            U0f=sum(U0fparts);
            
            %mean power frequency--don't use
            MPFreq(j,k)=U0f/U0;
            
        
        end
        %make resultant vector in transverse plane
        pathsxy(:,j)=sqrt((paths(:,1,j).^2)+(paths(:,2,j).^2)); 
        swaylengthxy(j)=swaylength(j,1)+swaylength(j,2);
        meandiffxy(j)=mean(abs(pathsxy(:,j)));
        rmserrorxy(j)=rms(pathsxy(:,j));
        swayarea_i=[];
        for m=1:(length(pathsxy)-1)
            swayarea_i(m)=abs((paths(m+1,2,j)*paths(m,1,j))-(paths(m,2,j)*paths(m+1,1,j))); %as per prieto 1996
        end
        swayarea(j)=(sum(swayarea_i))/(time*2);        
        %circle area(Prieto_1996) and ellipse area
        circlearea(j)=pi*(meandiffxy(j)+1.645*std(pathsxy(:,j)))^2;
        ellipsearea(j)=2*pi*3*sqrt(((std(paths(:,2,j)))^2*(std(paths(:,1,j)))^2)-((sum(paths(:,1,j).*paths(:,2,j)))/length(pathsxy))^2);
        meanfreqxy(j)=swaylengthxy(j)/(2*pi*time*meandiffxy(j));
        %Fractal dimension confidence circle
        FDCCxy(j)=(log10(length(pathsxy)))/(log10((length(pathsxy)*(2*(meandiffxy(j)+1.645*std(pathsxy(:,j)))))/swaylengthxy(j)));
        %Fractal dimension confidence ellipse
        FDCExy(j)=(log10(length(pathsxy)))/(log10((length(pathsxy)*sqrt((8*3*(sqrt(((std(paths(:,2,j)))^2*(std(paths(:,1,j)))^2)-((sum(paths(:,1,j).*paths(:,2,j)))/length(pathsxy))^2)))))/swaylengthxy(j)));
        %power spectral density computed using the sinusoidal multaper method with 8 tapers (as per prieto_1996)
%             nw=(8+1)/2; %8 tapers
        [Pdensxy(:,j),Pfreqxy(:,j)]=pmtm(pathsxy(:,j),nw,[],freq(j));
        %freq domain measures from 0.1465 to 5.01 Hz (Prieto does 0.15 to 5 Hz)
        %find indices closes to 0.15 hx and 5 Hx
        [MINlow MINlowi]=min(abs(Pfreqxy(:,j)-0.15));
        [MINhigh MINhighi]=min(abs(Pfreqxy(:,j)-5));
        TotPowerxy(j)=trapz(Pfreqxy(MINlowi:MINhighi,j),Pdensxy(MINlowi:MINhighi,j));
        %50 and 95% power frequency
        [min50 min50i]=min(abs(cumtrapz(Pfreqxy(MINlowi:MINhighi,j),Pdensxy(MINlowi:MINhighi,j))./TotPowerxy(j)-0.5));
        Powerxy50(j)=Pfreqxy(MINlowi+min50i);
        [min95 min95i]=min(abs(cumtrapz(Pfreqxy(MINlowi:MINhighi,j),Pdensxy(MINlowi:MINhighi,j))./TotPowerxy(j)-0.95));
        Powerxy95(j)=Pfreqxy(MINlowi+min95i);
        %keep in mind that units are in m, so values smaller than if were mm (multiply by 1000000--yes, if convert path to mm, then equivalent)
        %Tot Powerxy is redundant with RMS... although not if windowed from 0.15 hz to 5 hz like we did.
        
        %plot range or max point between distances
        Rangexy(j)=max(max(squareform(pdist([paths(:,1,j),paths(:,2,j)]))));
        
        %centroidal frequency and frequency dispersion as per Prieto, 1996 (equations 29 and 30)
        changeF=Pfreqxy(2,j)-Pfreqxy(1,j);
        hh=1;
        for h=MINlowi:MINhighi
            U2parts(hh)=(h*changeF)^2*Pdensxy(h,j);
            U1parts(hh)=(h*changeF)^1*Pdensxy(h,j);
            U0parts(hh)=Pdensxy(h,j);
            U0fparts(hh)=Pdens(h,j)*Pfreq(h,j);
            hh=hh+1;
        end
        U2=sum(U2parts);
        U1=sum(U1parts);
        U0=sum(U0parts);
        CFREQxy(j)=sqrt(U2/U0);
        FREQDxy(j)=sqrt(1-(U1^2)/(U2*U0));
        U0totalxy(j)=U0;
        U0f=sum(U0fparts);
        
        %U0 and totalpower measures should theoretically be the same.  They are not (probably because...
        %total power use trapz with respect to freq units.  Since RMS is analogous, no need to report.
        
        %Mean power frequency--don't use but good to have
        MPFreqxy(j)=U0f/U0;
        %Stochastic parameters, code written by V Lugade
        if convert=='y';
            try
                [SMP NSMP] = cop_smp_nsmp([paths(:,1,j),paths(:,2,j)], freq(j), 1); %1 in last input for graphing figures
                SMPdeltaTc(j)=SMP(1);
                SMPdeltaR2(j)=SMP(2);
                SMPDs(j)=SMP(3);
                SMPDl(j)=SMP(4);
                SMPHs(j)=SMP(5);
                SMPHl(j)=SMP(6);
                NSMPKs(j)=NSMP(1);
                NSMPKl(j)=NSMP(2);
                NSMPHs(j)=NSMP(3);
                NSMPHl(j)=NSMP(4);
            catch ME
                ME.message
                ME.stack.name
                ME.stack.line
                SMPdeltaTc(j)=0;
                SMPdeltaR2(j)=0;
                SMPDs(j)=0;
                SMPDl(j)=0;
                SMPHs(j)=0;
                SMPHl(j)=0;
                NSMPKs(j)=0;
                NSMPKl(j)=0;
                NSMPHs(j)=0;
                NSMPHl(j)=0;
                disp('stochastic parameter not calculated for...')
                if j==1
                    disp('...eyes open COP');
                elseif j==2
                    disp('...eyes open xCOP');
                elseif j==3
                    disp('...eyes open COM');
                elseif j==4
                    disp('...eyes open xCOM');
                elseif j==5
                    disp('...eyes closed COP');
                elseif j==6
                    disp('...eyes closed xCOP');
                elseif j==7
                    disp('...eyes closed COM');
                else j==8
                    disp('...eyes closed xCOM');
                end
                figure; %need to still have a figure for figure saving
            end
                    %Page Order, Eyes Open COP xCOP COM xCOM, then eyes closed COP xCOP COM xCOM, columns are x y
        else
            disp('stochastic parameters not calculated because units not in mm');
            SMPdeltaTc(j)=0;
            SMPdeltaR2(j)=0;
            SMPDs(j)=0;
            SMPDl(j)=0;
            SMPHs(j)=0;
            SMPHl(j)=0;
            NSMPKs(j)=0;
            NSMPKl(j)=0;
            NSMPHs(j)=0;
            NSMPHl(j)=0;
        end
    end
    %COP/COM sway ratio
    SwayRatioOPEN=swaylengthxy(1)/swaylengthxy(3);
    SwayRatioCLOSED=swaylengthxy(5)/swaylengthxy(7);
    
    %trouble shooting below
    %TotPowerCheck(i,1:2)=TotPower(1,1:2);
    %TotPowerxyCheck(i)=TotPowerxy(1);
    
    %Velocity of COP and COM variables
    for jj=1:4
        for k=1:2
            %root mean square velocity (RMS)
            rmserrorV(jj,k)=rms(vels(:,k,jj));
            %mean speed
            meanV(jj,k)=mean(abs(vels(:,k,jj)));
        end
        velsxy(:,jj)=sqrt((vels(:,1,jj).^2)+(vels(:,2,jj).^2));
        rmserrorVxy(jj)=rms(velsxy(:,jj));
        meanVxy(jj)=mean(abs(velsxy(:,jj)));
    end
    
    %set up for data export
    %no velocity variables (near end of list) for extrapolated variables
                SMPdeltaTc(j)=SMP(1);
            SMPdeltaR2(j)=SMP(2);
            SMPDs(j)=SMP(3);
            SMPDl(j)=SMP(4);
            SMPHs(j)=SMP(5);
            SMPHl(j)=SMP(6);
            NSMPKs(j)=NSMP(1);
            NSMPKl(j)=NSMP(2);
            NSMPHs(j)=NSMP(3);
            NSMPHl(j)=NSMP(4);
            
    V=1;
    newoutputOPEN_COP=[subject,1,V,AvgfootL,AvgfootW,height,Pheight,OPENfootalpha,OPENfootbeta,OPENbosL,...
        OPENbosW,OPENbosA,KickLimb,FirstUNIstance,MaxUNIboth,MaxUNIright,MaxUNIleft,MaxUNIkick,...
        MaxUNInonkick,OPENGRFasymm,OPENtrunksway(1),OPENtrunksway(2),OPENtrunkswayvel(1),...
        OPENtrunkswayvel(2),meanOPENpendlength,OPENPSDpercentauc(1),OPENPSDpercentauc(2),...
        OPENCOMposPSDpercentauc(1),OPENCOMposPSDpercentauc(2),OPENfreq,swaylengthxy(V),swaylength(V,1)...
        swaylength(V,2),meandiffxy(V),meandiff(V,1),meandiff(V,2),rmserrorxy(V),rmserror(V,1),...
        rmserror(V,2),Rangexy(V),maxdisp(V,1),maxdisp(V,2),meanfreqxy(V),meanfreq(V,1),meanfreq(V,2),...
        TotPowerxy(V),TotPower(V,1),TotPower(V,2),CFREQxy(V),CFREQ(V,1),CFREQ(V,2),FREQDxy(V),...
        FREQD(V,1),FREQD(V,2),U0totalxy(V),U0total(V,1),U0total(V,2),MPFreqxy(V),MPFreq(V,1),MPFreq(V,2),...
        swayarea(V),circlearea(V),ellipsearea(V),SwayRatioOPEN,rmserrorVxy(V),rmserrorV(V,1),rmserrorV(V,2),...
        meanVxy(V),meanV(V,1),meanV(V,2),SMPdeltaTc(V),SMPdeltaR2(V),SMPDs(V),SMPDl(V),SMPHs(V),...
        SMPHl(V),NSMPKs(V),NSMPKl(V),NSMPHs(V),NSMPHl(V)];
    
    V=2;
    newoutputOPEN_xCOP=[subject,1,V,AvgfootL,AvgfootW,height,Pheight,OPENfootalpha,OPENfootbeta,OPENbosL,...
        OPENbosW,OPENbosA,KickLimb,FirstUNIstance,MaxUNIboth,MaxUNIright,MaxUNIleft,MaxUNIkick,...
        MaxUNInonkick,OPENGRFasymm,OPENtrunksway(1),OPENtrunksway(2),OPENtrunkswayvel(1),...
        OPENtrunkswayvel(2),meanOPENpendlength,OPENPSDpercentauc(1),OPENPSDpercentauc(2),...
        OPENCOMposPSDpercentauc(1),OPENCOMposPSDpercentauc(2),OPENfreq,swaylengthxy(V),swaylength(V,1)...
        swaylength(V,2),meandiffxy(V),meandiff(V,1),meandiff(V,2),rmserrorxy(V),rmserror(V,1),...
        rmserror(V,2),Rangexy(V),maxdisp(V,1),maxdisp(V,2),meanfreqxy(V),meanfreq(V,1),meanfreq(V,2),...
        TotPowerxy(V),TotPower(V,1),TotPower(V,2),CFREQxy(V),CFREQ(V,1),CFREQ(V,2),FREQDxy(V),...
        FREQD(V,1),FREQD(V,2),U0totalxy(V),U0total(V,1),U0total(V,2),MPFreqxy(V),MPFreq(V,1),MPFreq(V,2),...
        swayarea(V),circlearea(V),ellipsearea(V),SwayRatioOPEN,0,0,0,0,0,0,SMPdeltaTc(V),SMPdeltaR2(V),SMPDs(V),SMPDl(V),SMPHs(V),...
        SMPHl(V),NSMPKs(V),NSMPKl(V),NSMPHs(V),NSMPHl(V)];
    
    V=3;
    newoutputOPEN_COM=[subject,1,V,AvgfootL,AvgfootW,height,Pheight,OPENfootalpha,OPENfootbeta,OPENbosL,...
        OPENbosW,OPENbosA,KickLimb,FirstUNIstance,MaxUNIboth,MaxUNIright,MaxUNIleft,MaxUNIkick,...
        MaxUNInonkick,OPENGRFasymm,OPENtrunksway(1),OPENtrunksway(2),OPENtrunkswayvel(1),...
        OPENtrunkswayvel(2),meanOPENpendlength,OPENPSDpercentauc(1),OPENPSDpercentauc(2),...
        OPENCOMposPSDpercentauc(1),OPENCOMposPSDpercentauc(2),OPENfreq,swaylengthxy(V),swaylength(V,1)...
        swaylength(V,2),meandiffxy(V),meandiff(V,1),meandiff(V,2),rmserrorxy(V),rmserror(V,1),...
        rmserror(V,2),Rangexy(V),maxdisp(V,1),maxdisp(V,2),meanfreqxy(V),meanfreq(V,1),meanfreq(V,2),...
        TotPowerxy(V),TotPower(V,1),TotPower(V,2),CFREQxy(V),CFREQ(V,1),CFREQ(V,2),FREQDxy(V),...
        FREQD(V,1),FREQD(V,2),U0totalxy(V),U0total(V,1),U0total(V,2),MPFreqxy(V),MPFreq(V,1),MPFreq(V,2),...
        swayarea(V),circlearea(V),ellipsearea(V),SwayRatioOPEN,rmserrorVxy(2),rmserrorV(2,1),rmserrorV(2,2),...
        meanVxy(2),meanV(2,1),meanV(2,2),SMPdeltaTc(V),SMPdeltaR2(V),SMPDs(V),SMPDl(V),SMPHs(V),...
        SMPHl(V),NSMPKs(V),NSMPKl(V),NSMPHs(V),NSMPHl(V)];
    
    V=4;
    newoutputOPEN_xCOM=[subject,1,V,AvgfootL,AvgfootW,height,Pheight,OPENfootalpha,OPENfootbeta,OPENbosL,...
        OPENbosW,OPENbosA,KickLimb,FirstUNIstance,MaxUNIboth,MaxUNIright,MaxUNIleft,MaxUNIkick,...
        MaxUNInonkick,OPENGRFasymm,OPENtrunksway(1),OPENtrunksway(2),OPENtrunkswayvel(1),...
        OPENtrunkswayvel(2),meanOPENpendlength,OPENPSDpercentauc(1),OPENPSDpercentauc(2),...
        OPENCOMposPSDpercentauc(1),OPENCOMposPSDpercentauc(2),OPENfreq,swaylengthxy(V),swaylength(V,1)...
        swaylength(V,2),meandiffxy(V),meandiff(V,1),meandiff(V,2),rmserrorxy(V),rmserror(V,1),...
        rmserror(V,2),Rangexy(V),maxdisp(V,1),maxdisp(V,2),meanfreqxy(V),meanfreq(V,1),meanfreq(V,2),...
        TotPowerxy(V),TotPower(V,1),TotPower(V,2),CFREQxy(V),CFREQ(V,1),CFREQ(V,2),FREQDxy(V),...
        FREQD(V,1),FREQD(V,2),U0totalxy(V),U0total(V,1),U0total(V,2),MPFreqxy(V),MPFreq(V,1),MPFreq(V,2),...
        swayarea(V),circlearea(V),ellipsearea(V),SwayRatioOPEN,0,0,0,0,0,0,SMPdeltaTc(V),SMPdeltaR2(V),SMPDs(V),SMPDl(V),SMPHs(V),...
        SMPHl(V),NSMPKs(V),NSMPKl(V),NSMPHs(V),NSMPHl(V)];
    
    V=5;
    newoutputCLOSED_COP=[subject,2,V,AvgfootL,AvgfootW,height,Pheight,CLOSEDfootalpha,CLOSEDfootbeta,CLOSEDbosL,...
        CLOSEDbosW,CLOSEDbosA,KickLimb,FirstUNIstance,MaxUNIboth,MaxUNIright,MaxUNIleft,MaxUNIkick,...
        MaxUNInonkick,CLOSEDGRFasymm,CLOSEDtrunksway(1),CLOSEDtrunksway(2),CLOSEDtrunkswayvel(1),...
        CLOSEDtrunkswayvel(2),meanCLOSEDpendlength,CLOSEDPSDpercentauc(1),CLOSEDPSDpercentauc(2),...
        CLOSEDCOMposPSDpercentauc(1),CLOSEDCOMposPSDpercentauc(2),CLOSEDfreq,swaylengthxy(V),swaylength(V,1)...
        swaylength(V,2),meandiffxy(V),meandiff(V,1),meandiff(V,2),rmserrorxy(V),rmserror(V,1),...
        rmserror(V,2),Rangexy(V),maxdisp(V,1),maxdisp(V,2),meanfreqxy(V),meanfreq(V,1),meanfreq(V,2),...
        TotPowerxy(V),TotPower(V,1),TotPower(V,2),CFREQxy(V),CFREQ(V,1),CFREQ(V,2),FREQDxy(V),...
        FREQD(V,1),FREQD(V,2),U0totalxy(V),U0total(V,1),U0total(V,2),MPFreqxy(V),MPFreq(V,1),MPFreq(V,2),...
        swayarea(V),circlearea(V),ellipsearea(V),SwayRatioCLOSED,rmserrorVxy(3),rmserrorV(3,1),rmserrorV(3,2),...
        meanVxy(3),meanV(3,1),meanV(3,2),SMPdeltaTc(V),SMPdeltaR2(V),SMPDs(V),SMPDl(V),SMPHs(V),...
        SMPHl(V),NSMPKs(V),NSMPKl(V),NSMPHs(V),NSMPHl(V)];
    
    V=6;
    newoutputCLOSED_xCOP=[subject,2,V,AvgfootL,AvgfootW,height,Pheight,CLOSEDfootalpha,CLOSEDfootbeta,CLOSEDbosL,...
        CLOSEDbosW,CLOSEDbosA,KickLimb,FirstUNIstance,MaxUNIboth,MaxUNIright,MaxUNIleft,MaxUNIkick,...
        MaxUNInonkick,CLOSEDGRFasymm,CLOSEDtrunksway(1),CLOSEDtrunksway(2),CLOSEDtrunkswayvel(1),...
        CLOSEDtrunkswayvel(2),meanCLOSEDpendlength,CLOSEDPSDpercentauc(1),CLOSEDPSDpercentauc(2),...
        CLOSEDCOMposPSDpercentauc(1),CLOSEDCOMposPSDpercentauc(2),CLOSEDfreq,swaylengthxy(V),swaylength(V,1)...
        swaylength(V,2),meandiffxy(V),meandiff(V,1),meandiff(V,2),rmserrorxy(V),rmserror(V,1),...
        rmserror(V,2),Rangexy(V),maxdisp(V,1),maxdisp(V,2),meanfreqxy(V),meanfreq(V,1),meanfreq(V,2),...
        TotPowerxy(V),TotPower(V,1),TotPower(V,2),CFREQxy(V),CFREQ(V,1),CFREQ(V,2),FREQDxy(V),...
        FREQD(V,1),FREQD(V,2),U0totalxy(V),U0total(V,1),U0total(V,2),MPFreqxy(V),MPFreq(V,1),MPFreq(V,2),...
        swayarea(V),circlearea(V),ellipsearea(V),SwayRatioCLOSED,0,0,0,0,0,0,SMPdeltaTc(V),SMPdeltaR2(V),SMPDs(V),SMPDl(V),SMPHs(V),...
        SMPHl(V),NSMPKs(V),NSMPKl(V),NSMPHs(V),NSMPHl(V)];
    
    V=7;
    newoutputCLOSED_COM=[subject,2,V,AvgfootL,AvgfootW,height,Pheight,CLOSEDfootalpha,CLOSEDfootbeta,CLOSEDbosL,...
        CLOSEDbosW,CLOSEDbosA,KickLimb,FirstUNIstance,MaxUNIboth,MaxUNIright,MaxUNIleft,MaxUNIkick,...
        MaxUNInonkick,CLOSEDGRFasymm,CLOSEDtrunksway(1),CLOSEDtrunksway(2),CLOSEDtrunkswayvel(1),...
        CLOSEDtrunkswayvel(2),meanCLOSEDpendlength,CLOSEDPSDpercentauc(1),CLOSEDPSDpercentauc(2),...
        CLOSEDCOMposPSDpercentauc(1),CLOSEDCOMposPSDpercentauc(2),CLOSEDfreq,swaylengthxy(V),swaylength(V,1)...
        swaylength(V,2),meandiffxy(V),meandiff(V,1),meandiff(V,2),rmserrorxy(V),rmserror(V,1),...
        rmserror(V,2),Rangexy(V),maxdisp(V,1),maxdisp(V,2),meanfreqxy(V),meanfreq(V,1),meanfreq(V,2),...
        TotPowerxy(V),TotPower(V,1),TotPower(V,2),CFREQxy(V),CFREQ(V,1),CFREQ(V,2),FREQDxy(V),...
        FREQD(V,1),FREQD(V,2),U0totalxy(V),U0total(V,1),U0total(V,2),MPFreqxy(V),MPFreq(V,1),MPFreq(V,2),...
        swayarea(V),circlearea(V),ellipsearea(V),SwayRatioCLOSED,rmserrorVxy(4),rmserrorV(4,1),rmserrorV(4,2),...
        meanVxy(4),meanV(4,1),meanV(4,2),SMPdeltaTc(V),SMPdeltaR2(V),SMPDs(V),SMPDl(V),SMPHs(V),...
        SMPHl(V),NSMPKs(V),NSMPKl(V),NSMPHs(V),NSMPHl(V)];
    
    V=8;
    newoutputCLOSED_xCOM=[subject,2,V,AvgfootL,AvgfootW,height,Pheight,CLOSEDfootalpha,CLOSEDfootbeta,CLOSEDbosL,...
        CLOSEDbosW,CLOSEDbosA,KickLimb,FirstUNIstance,MaxUNIboth,MaxUNIright,MaxUNIleft,MaxUNIkick,...
        MaxUNInonkick,CLOSEDGRFasymm,CLOSEDtrunksway(1),CLOSEDtrunksway(2),CLOSEDtrunkswayvel(1),...
        CLOSEDtrunkswayvel(2),meanCLOSEDpendlength,CLOSEDPSDpercentauc(1),CLOSEDPSDpercentauc(2),...
        CLOSEDCOMposPSDpercentauc(1),CLOSEDCOMposPSDpercentauc(2),CLOSEDfreq,swaylengthxy(V),swaylength(V,1)...
        swaylength(V,2),meandiffxy(V),meandiff(V,1),meandiff(V,2),rmserrorxy(V),rmserror(V,1),...
        rmserror(V,2),Rangexy(V),maxdisp(V,1),maxdisp(V,2),meanfreqxy(V),meanfreq(V,1),meanfreq(V,2),...
        TotPowerxy(V),TotPower(V,1),TotPower(V,2),CFREQxy(V),CFREQ(V,1),CFREQ(V,2),FREQDxy(V),...
        FREQD(V,1),FREQD(V,2),U0totalxy(V),U0total(V,1),U0total(V,2),MPFreqxy(V),MPFreq(V,1),MPFreq(V,2),...
        swayarea(V),circlearea(V),ellipsearea(V),SwayRatioCLOSED,0,0,0,0,0,0,SMPdeltaTc(V),SMPdeltaR2(V),SMPDs(V),SMPDl(V),SMPHs(V),...
        SMPHl(V),NSMPKs(V),NSMPKl(V),NSMPHs(V),NSMPHl(V)];
    
    
    newoutputA=[newoutputOPEN_COP; newoutputOPEN_xCOP; newoutputOPEN_COM; newoutputOPEN_xCOM; ...
        newoutputCLOSED_COP; newoutputCLOSED_xCOP; newoutputCLOSED_COM; newoutputCLOSED_xCOM];
    
    newoutput=[newoutput;newoutputA];
    
    %Save figures
    if fileoutput=='y'
        figure(1)
        saveas(gcf,[outfolder,'figures\','FAFRa',num2str(i,'%03d'),'eyesopenBOS.fig'],'fig');
        figure(2)
        saveas(gcf,[outfolder,'figures\','FAFRa',num2str(i,'%03d'),'eyesclosedBOS.fig'],'fig');
        figure(3)
        saveas(gcf,[outfolder,'figures\','FAFRa',num2str(i,'%03d'),'eyesopen.fig'],'fig');
        figure(4)
        saveas(gcf,[outfolder,'figures\','FAFRa',num2str(i,'%03d'),'eyesclosed.fig'],'fig');
        figure(5)
        saveas(gcf,[outfolder,'figures\','FAFRa',num2str(i,'%03d'),'stochasticCOPeyesopen.fig'],'fig');
        figure(6)
        saveas(gcf,[outfolder,'figures\','FAFRa',num2str(i,'%03d'),'stochasticxCOPeyesopen.fig'],'fig');
        figure(7)
        saveas(gcf,[outfolder,'figures\','FAFRa',num2str(i,'%03d'),'stochasticCOMeyesopen.fig'],'fig');
        figure(8)
        saveas(gcf,[outfolder,'figures\','FAFRa',num2str(i,'%03d'),'stochasticxCOMeyesopen.fig'],'fig');
        figure(9)
        saveas(gcf,[outfolder,'figures\','FAFRa',num2str(i,'%03d'),'stochasticCOPeyesclosed.fig'],'fig');
        figure(10)
        saveas(gcf,[outfolder,'figures\','FAFRa',num2str(i,'%03d'),'stochasticxCOPeyesclosed.fig'],'fig');
        figure(11)
        saveas(gcf,[outfolder,'figures\','FAFRa',num2str(i,'%03d'),'stochasticCOMeyesclosed.fig'],'fig');
        figure(12)
        saveas(gcf,[outfolder,'figures\','FAFRa',num2str(i,'%03d'),'stochasticxCOMeyesclosed.fig'],'fig');
    end
end %end loop for each subject

% %for first line of spreadsheet
% header={'subject','eyes','variable','AvgfootL','AvgfootW','height','Pheight','footalpha','footbeta',...
%     'bosL','bosW','bosA','KickLimb','FirstUNIstance','MaxUNIboth','MaxUNIright','MaxUNIleft','MaxUNIkick',...
%     'MaxUNInonkick','GRFasymm','trunkswayX','trunkswayY','trunkswayvelX','trunkswayvelY','meanpendlength',...
%     'COPPSDpercentaucX','COPPSDpercentaucY','COMPSDpercentaucX','COMPSDpercentaucY','freq','SwayLengthXY',...
%     'SwayLengthX','SwayLengthY','meandiffXY','meandiffX','meandiffY','rmserrorXY','rmserrorX','rmserrorY',...
%     'Rangexy','maxdispX','maxdispY','meanfreqXY','meanfreqX','meanfreqY','TotPowerXY','TotPowerX','TotPowerY',...
%     'CFREQXY','CFREQX','CFREQY','FREQDXY','FREQDX','FREQDY','U0totalXY','U0totalX','U0totalY',...
%     'MPFreqXY','MPFreqX','MPFreqY','swayarea','circlearea','ellipsearea','SwayRatio',...
%     'rmserrorVXY','rmserrorVX','rmserrorVY','meanVXY','meanVX','meanVY','SMPdeltaTc','SMPdeltaR2','SMPDs'...
%       'SMPDl','SMPHs','SMPHl','NSMPKs','NSMPK1','NSMPHs','NSMPHl'}



%Identify if user selected existing file or created new file
fileexist=exist(outfile);
if fileexist==2 %file does exist
    [N,T,output2]=xlsread(outfile);
    output=[N;newoutput];
else
    output=newoutput;
end

% output={header;output};
outputsize=size(output);
range=['A2:CB',num2str(outputsize(1)+1)];

if fileoutput=='y'
    xlswrite(outfile,output,range);
end




%Export and append to 
%I:\INVIVO\RESEARCH\FAFR_11-006984_Amin\MOTIONVIDEO\CompiledData\Sway
%FAFR_bilateralsway_compiled.xls
%use outputfile boolean in beginning for option to not update forceplate

%Export images of COPnet xCOP and COM and xCOM global

%clear variables before each i iteration.

%Variables:

%stochastic parameters (Chiari 2002 gives additanl references)
%short-term scaling exponent (Chiari 2002)
%short-term diffusion coefficient (Chiari 2002, Maurer 2005)
%long-term scaling exponent (Chiari 2002)
%long-term diffusion coefficent (Chiari 2002, Maurer 2005)
%diffusion coefficient
%time lag corresponding to a pure random motion (Chiari 2002)
%Hurst rescaled range analysis (Lin 2008)
%detrended fluction analysis (Lin 2008)

%critical point coordinate from the sdf time (, Maurer 2005)
%critical point coordinate from the sdf amplitude (, Maurer 2005)

%spectral entropy, approximate entropy, and singular value decomposition
%spectrum entropy (Sabatini 2000)

%Tests relative to BOS length, width, area?
%Also look at ratio of eyes closed / eyes open (Romber ratio Sabatini 2000)
