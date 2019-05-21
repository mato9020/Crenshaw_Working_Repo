function [COPdist,COPds,COPfilt,VAR] = COPSwayCheck(endTime,fileName)
%UNTITLED Summary of this function goes here
%   AP- Anterior Posterior (x)
%   ML- Medial Lateral (y)
%   TV- transverse (xy)
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 

int16 endTime;
char fileName; 
uint8 SRate; %Sampling rate of Force plate determined by trial
swayAreaM = [];
% [AP,ML] = COPplanes;



%Below lines are for CORTEX based .forces files  % commented out
% INPUT=importdata(fileName,'\t');
% COP=INPUT.data(:,5:6);
%Below lines are for QTM TSV files
INPUT=importdata(fileName,'\t',27);
COP=INPUT.data(:,7:8);

%Sampling rate adjustments to match lowest sample rate
COP=COP(1:endTime*1200,:); %assumes 1200 Hz, then we only want 30 s
COPds=downsample(COP,2); %downsamples to 120 Hz if set to "10".  downsamples to 600 Hz if set to 2.
COfreq=10; %Cut-off freq for lowpass filter of COP and GRF




[z,p] = butter(2,COfreq/(600/2));  %Assumes sampling frequency is 600 Hz, 2nd order, or forth order recursive
COPfilt(:,1)=filtfilt(z,p,COPds(:,1));
COPfilt(:,2)=filtfilt(z,p,COPds(:,2));

AP=(COPfilt(:,1)-mean(COPfilt(:,1)));ML=(COPfilt(:,2)-mean(COPfilt(:,2)));%differentiates each COP signal and assigns x and y to AP and ML arrays respectively
COPdist = hypot((AP),(ML)); %RD time series (Prietto)
stdAP = std(AP); stdML = std(ML); %standard deviation in each plane


%Variable outputs
APrms = rms(AP);MLrms = rms(ML);
TOTEX = sum(diff(COPdist)); %total path length
TOTEXap = sum(diff(AP));
TOTEXml = sum(diff(ML));
% MDIST = mean(COPdist); %Mean Transverse distance
MDIST = mean(COPdist);
RDIST = rms(COPdist);
MVELO = TOTEX/endTime; MVELOap = TOTEXap/endTime; MVELOml = TOTEXml/endTime;
stdTV = (RDIST^2-MDIST^2)^(1/2);
ccrad = MDIST+1.645*std(COPdist);
CCArea = pi*(ccrad)^2;

CEArea = 6*pi*((stdAP^2).*(stdML^2)-(sum(AP.*ML)/length(ML))^2)^(1/2);

for m=1:length(AP)-1
    swayAreaM(m) = abs(AP(m+1)*ML(m)-AP(m)*ML(m+1));
end
swayArea = 1/(2*endTime)*(sum(swayAreaM));

MFREQ = MVELO/(2*pi*MDIST);
VAR={'MDIST','stdTV','stdAP','stdML','Transverse rms','APrms','MLrms','TOTEX','MVELO','CCArea','CEArea','swayArea','MFREQ','ccrad';MDIST,stdTV,stdAP,stdML,RDIST,APrms,MLrms,TOTEX,MVELO,CCArea,CEArea,swayArea,MFREQ,ccrad};


plot(COPds(:,1)-mean(COPds(:,1)),COPds(:,2)-mean(COPds(:,2)))
hold on
plot(AP,ML,'r')

end

