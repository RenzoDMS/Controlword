close all
clear all
clc


%addpath('C:\Users\lapa-8\Documents\Philemon\MATLAB')%MODIFY according to
%the folder where the files to be read, the script and the functions are.
[filename, pathname_ascii] = uigetfile({'*.txt'},'Select classic ADCP ASCII output from WRII');
fullName=[pathname_ascii filename]; %Enter here name of classic ASCII output file
%from WinRiver you wish to read
% fullName = 'Maca000_ASC.TXT';
ignoreBS=0;
screenData=1;
A=readADCPfile(fullName,screenData,ignoreBS);

%% Display warning
% wrnBox= warndlg('A transect with lots of missing values/bad bins will yield erronous flux calculations due to uncertainties in the extrapolation');
% uiwait(wrnBox);


% Get name and date of station
%% Enter name of station and date of measure.
stationname = inputdlg('Please enter the name of the station:', 'Station');
stationname = num2str(stationname{1});
shore =  questdlg('Please define the starting bank of the transect used:', 'Starting shore', 'left', 'right','left');
shore = [shore ' bank [m]'];
yr = num2str(A.Sup.year(1));
month = num2str(A.Sup.month(1));
day = num2str(A.Sup.day(1));
date = [day,'-', month,'-', yr];

namedate = [stationname, ' on the ', date];
%% Extract ADCP data required
% Get date of measure:
yr = A.Sup.year(1);
month = A.Sup.month(1);
day = A.Sup.day(1);
date = [day month yr];
date = num2str(date);
namedate = [stationname, ' on the ', date];
%EnsN=A.Sup.ensNo;
% Beam intensity in counts
Intensity1=A.Wat.backscatter; %This is a 3D matrice with each slice corresponding to BS of one beam
IntB1 = Intensity1(:,:,1);% Beam intensity matrices
IntB2 = Intensity1(:,:,2);
IntB3 = Intensity1(:,:,3);
IntB4 = Intensity1(:,:,4);
indb1 = find(IntB1 == 255); %Not measured VALUE of intensity (COUNTS)
IntB1(indb1)=NaN;
indb2 = find(IntB2 == 255); %Not measured VALUE of intensity (COUNTS)
IntB2(indb2)=NaN;
indb3 = find(IntB3 == 255); %Not measured VALUE of intensity (COUNTS)
IntB3(indb3)=NaN;
indb4 = find(IntB4 == 255); %Not measured VALUE of intensity (COUNTS)
IntB4(indb4)=NaN;
IntBav = (IntB1 + IntB2 +IntB3 +IntB4)./4;
IntB1B2 = (IntB1 +IntB2)./2;
IntB1B3 = (IntB1 +IntB3)./2;
IntB1B4 = (IntB1 +IntB4)./2;
IntB2B3 = (IntB2+IntB3)./2;
IntB2B4 =(IntB2+IntB4)./2;
IntB3B4 = (IntB3+IntB4)./2;

% Velocity
vel_dir = A.Wat.vDir;     %velocity direction
velE = A.Wat.vEast./100;  %east velocity
velN = A.Wat.vNorth./100; %north velocity
velup = A.Wat.vVert./100; %vertical velocity
velMag = A.Wat.vMag./100; %velocity magnitude

%Range vector to give the depth of each bin
Rangesvec = A.Wat.binDepth(:,1);


DistE = A.Nav.totDistEast;
DistN = A.Nav.totDistNorth;
% Beam depths
% Test to see if correction of 20 deg angle works better

% beam_angle=20;
% DepthB1 = A.Nav.depth(:,1)./cos(beam_angle/180*pi);
% DepthB2 = A.Nav.depth(:,2)./cos(beam_angle/180*pi);
% DepthB3 = A.Nav.depth(:,3)./cos(beam_angle/180*pi);
% DepthB4 = A.Nav.depth(:,4)./cos(beam_angle/180*pi);
DepthB1 = A.Nav.depth(:,1);
DepthB2 = A.Nav.depth(:,2);
DepthB3 = A.Nav.depth(:,3);
DepthB4 = A.Nav.depth(:,4);
DepthBave = (DepthB1+DepthB2+DepthB3+DepthB4)./4;
DepthB1B2 = (DepthB1 +DepthB2)./2;
DepthB1B3 = (DepthB1 +DepthB3)./2;
DepthB1B4 = (DepthB1 +DepthB4)./2;
DepthB2B3 = (DepthB2+DepthB3)./2;
DepthB2B4 =(DepthB2+DepthB4)./2;
DepthB3B4 = (DepthB3+DepthB4)./2;

% Total distance made good
DMG = A.Nav.dmg;

%Temperature
Temp = A.Sensor.temp_degC;

% Number of bins
Nbins=A.Sup.bins;
EnsN=A.Sup.ensNo;
Nens=size(EnsN,1);
indStart=1;
indEnd=Nens;
%% Added SCRiPT SPECIFIC VALUES
%Number of cells to add
% BinSize = double(A.Sup.binSize_cm)/100; % CONVERT VALUE FROM UINT32 TO DOUBLE

BinSize = double(A.Sup.binSize_cm)/100; % CONVERT VALUE FROM UINT32 TO DOUBLE
cells2add = fix(A.Wat.binDepth(1,1))/BinSize;
MatrixUpperPart = nan(cells2add,A.Sup.noe);

noe=A.Sup.noe;
Depth_bins = A.Wat.binDepth;
Q_units = A.Q.unit;


MaxDepth=max(DepthBave);
[nligtemp, ncol] = size(Depth_bins);
% Define depths of cells added to reach the water free surface
TopBinDepths=MatrixUpperPart;
for i = 1: cells2add
    TopBinDepths(i,:) =  Depth_bins(1,1)- double((cells2add+1-i).*BinSize);
end

%Add cells at the bottom
if MaxDepth>max(Depth_bins(:,1))
    cells2add_bot = ceil(MaxDepth - max(Depth_bins(:,1)))/BinSize+1; % we add 1 for security... doesnt matter as converted to nan later
    MatrixBotPart = nan(cells2add_bot, A.Sup.noe);
    BotBinDepths = MatrixBotPart;
    for i = 1: cells2add_bot
        BotBinDepths(i,:) =  Depth_bins(end,1)+ (i*BinSize);
    end
else BotBinDepths = [];
    MatrixBotPart =[];
end

Depth_binsTot= [TopBinDepths;Depth_bins;BotBinDepths];

Depth_Rouse = Depth_binsTot;
%% Treat unreal Q values
% Remove not real flow values (<0 or unmeasured)
ind_notmeasQ = find(Q_units == 2147483647);
ind_notRealQ= find(Q_units<=0);
Q_units(ind_notRealQ)=NaN;
Q_units(ind_notmeasQ)=NaN;
Q2recons= [MatrixUpperPart;Q_units;MatrixBotPart; ];
vel2recons= [MatrixUpperPart;velMag;MatrixBotPart ];


%% Find indices which correspond to max Depth for each vertical
%ind_maxDepth=nan(2,noe);
%% INTERPOLATE DEPTHBAVE TO HAVE NO NAN
indepthBave=find(DepthBave<=0);
DepthBave(indepthBave)=NaN;
for i = 2:numel(DepthBave)
    if abs(DepthBave(i)-DepthBave(i-1))> 1
        DepthBave(i) = DepthBave(i-1);
    end
end
DepthBave = inpaint_nans(DepthBave,2);
%% --------------------
% Find indices of last measured value
nlin_depth = nan(size(DepthBave,1),1);
ncol_depth = nan(size(DepthBave,1),1);

for i =1:noe
    [nlin_depth(i), ncol_depth(i)]= find(Depth_binsTot(:,i)<=DepthBave(i) + BinSize & Depth_binsTot(:,i)>=DepthBave(i)-BinSize, 1,'Last'); %N.B. A.Q.endDepth is the last measured value
end
%for i =1:noe
    %[nlin_depth(i), ncol_depth(i)]= find(Depth_binsTot(:,i)<=DepthBave(i) + 0.5 & Depth_binsTot(:,i)>=DepthBave(i)-0.5, 1,'Last'); %N.B. A.Q.endDepth is the last measured value
%end

for i = 1:noe
    Depth_Rouse(nlin_depth(i):end,i)=NaN; % Convert back unmeasureable values to Nan (below the bed)
    %The added (DepthBave(i)-endDepth(i))/BinSize corresponds to the
    %difference in bin number between last measurement and actual bed depth
end
% Arrange bed profile
DepthBave_Plot = DepthBave;
indDepth= find(DepthBave_Plot<=0); %Remove none measured values to
% have a smoother bed profile
DepthBave_Plot(indDepth)=NaN;
%% Try to remove peaks from DepthBav (smoothen the bed)
for i = 2:numel(DepthBave_Plot)
    if abs(DepthBave_Plot(i)-DepthBave_Plot(i-1))> 1
        DepthBave_Plot(i) = DepthBave_Plot(i-1);
    end
end

DepthBave_plot = inpaint_nans(DepthBave_Plot,2);
endDepth_Plot = A.Q.endDepth;
indendDepth= find(endDepth_Plot<=0); %Remove none measured values to
%% Plot average intensity to see number of verticals

% have a smoother bed profile
endDepth_Plot(indendDepth)=NaN;
figure;
pcolor(DMG,Rangesvec,IntB1);
hold on
shading interp
colormap(jet)
t=colorbar;
xlabel('Distance [m]');
set(gca,'YDir','reverse');
set(get(t,'ylabel'),'String', 'Backscatter [dB]');
ylabel('Depth [m]');
title(sprintf('Measured backscatter at %s',...
    namedate),'FontSize', 10, 'FontWeight','bold')
hold on;
plot(DMG, DepthBave_Plot, 'color', 'black', 'LineWidth',1.5)
hold on;
plot(DMG,endDepth_Plot, 'color', 'black', 'LineWidth',0.75)
hold off;

%% -----------------------------------------------
% CONC CALC USING ROUSE PROFILE
% --------------------------------------------
Conc_coarse = nan(size(Depth_Rouse,1), size(Depth_Rouse,2));
Conc_fine = nan(size(Depth_Rouse,1), size(Depth_Rouse,2));

inputOptions={' One  ', ' Two ', 'Three', 'Four','Five'};
bttnsOredring = [1,5];
defSelection=inputOptions{1};
numVert=bttnChoiceDialog(inputOptions, 'Verticals for concentrations', defSelection,...
    'Select the number of verticals you have to define the concentration of fines:',bttnsOredring);
% fprintf( 'User selection "%s"\n', inputOptions{numVert});
switch numVert
    case 1
        %% CASE WITH ONE VERTICAL.
        promptV1={sprintf('Distance from %s', shore),'Enter WASHLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter WASHLOAD Rouse number (W):', 'Enter BEDLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter BEDLOAD Rouse number (W):'};
        name_promptV1='Vertical 1';
        numlines=1;
        onevert=inputdlg(promptV1,name_promptV1,numlines);
        distV1=str2double(onevert(1)); %V1 Distance
        Cz0V1=str2double(onevert(2)); %Base concentration  WL
        W_V1 = str2double(onevert(3)); %Rouse exposant WL
        Cz02V1=str2double(onevert(4)); %Base concentration BL
        W2_V1 = str2double(onevert(5)); %Rouse exposant BL
        
        indV1 = find(abs((round(DMG(:))- distV1))<= 1.5,1, 'first' );
        for i = 1:noe %Loop on verticals
            h = max(Depth_Rouse(:,i))+0.5;
            z0 = h-Depth_Rouse(:,i);
            a = 0.05*h;
            for j = 1:numel(Depth_Rouse(:,1)) %Loop on bins
                Conc_fine(j,i) =  Cz0V1.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W_V1;
                Conc_coarse(j,i) = Cz02V1.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W2_V1;
            end
        end
    case 2
        %% CASE WITH TWO VERTICALS.
        promptV1={sprintf('Distance from %s', shore),'Enter WASHLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter WASHLOAD Rouse number (W):', 'Enter BEDLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter BEDLOAD Rouse number (W):'};
        name_promptV1='Vertical 1';
        numlines=1;
        onevert=inputdlg(promptV1,name_promptV1,numlines);
        distV1=str2double(onevert(1)); %V1 Distance
        Cz0V1=str2double(onevert(2)); %Base concentration  WL
        W_V1 = str2double(onevert(3)); %Rouse exposant WL
        Cz02V1=str2double(onevert(4)); %Base concentration BL
        W2_V1 = str2double(onevert(5)); %Rouse exposant BL
        %2nd vertical
        promptV2={sprintf('Distance from %s', shore),'Enter WASHLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter WASHLOAD Rouse number (W):', 'Enter BEDLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter BEDLOAD Rouse number (W):'};
        name_promptV2='Vertical 2';
        numlines=1;
        onevert=inputdlg(promptV2,name_promptV2,numlines);
        distV2=str2double(onevert(1)); %V1 Distance
        Cz0V2=str2double(onevert(2)); %Base concentration  WL
        W_V2 = str2double(onevert(3)); %Rouse exposant WL
        Cz02V2=str2double(onevert(4)); %Base concentration BL
        W2_V2 = str2double(onevert(5)); %Rouse exposant BL
        
        indV1 = find(abs((round(DMG(:))- distV1))<= 1.5,1, 'first' );
        indV2 = find(abs((round(DMG(:))- distV2))<= 1.5,1, 'first' );
        
        for i = 1:floor((indV1+indV2)/2) %Loop on verticals
            h = max(Depth_Rouse(:,i))+0.5;
            z0 = h-Depth_Rouse(:,i);
            a = 0.05*h;
            for j = 1:numel(Depth_Rouse(:,1)) %Loop on bins
                Conc_fine(j,i) =  Cz0V1.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W_V1;
                Conc_coarse(j,i) = Cz02V1.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W2_V1;
            end
            
        end
        
        for i = floor((indV1+indV2)/2)+1:noe %Loop on verticals
            h = max(Depth_Rouse(:,i))+0.5;
            z0 = h-Depth_Rouse(:,i);
            a = 0.05*h;
            for j = 1:numel(Depth_Rouse(:,1)) %Loop on bins
                Conc_fine(j,i) =  Cz0V2.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W_V2;
                Conc_coarse(j,i) = Cz02V2.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W2_V2;
            end
        end
    case 3
        %% CASE WITH THREE VERTICALS.
        promptV1={sprintf('Distance from %s', shore),'Enter WASHLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter WASHLOAD Rouse number (W):', 'Enter BEDLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter BEDLOAD Rouse number (W):'};
        name_promptV1='Vertical 1';
        numlines=1;
        onevert=inputdlg(promptV1,name_promptV1,numlines);
        distV1=str2double(onevert(1)); %V1 Distance
        Cz0V1=str2double(onevert(2)); %Base concentration  WL
        W_V1 = str2double(onevert(3)); %Rouse exposant WL
        Cz02V1=str2double(onevert(4)); %Base concentration BL
        W2_V1 = str2double(onevert(5)); %Rouse exposant BL
        %2nd vertical
        promptV2={sprintf('Distance from %s', shore),'Enter WASHLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter WASHLOAD Rouse number (W):', 'Enter BEDLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter BEDLOAD Rouse number (W):'};
        name_promptV2='Vertical 2';
        numlines=1;
        onevert=inputdlg(promptV2,name_promptV2,numlines);
        distV2=str2double(onevert(1)); %V1 Distance
        Cz0V2=str2double(onevert(2)); %Base concentration  WL
        W_V2 = str2double(onevert(3)); %Rouse exposant WL
        Cz02V2=str2double(onevert(4)); %Base concentration BL
        W2_V2 = str2double(onevert(5)); %Rouse exposant BL
        % 3rd vertical
        promptV3={sprintf('Distance from %s', shore),'Enter WASHLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter WASHLOAD Rouse number (W):', 'Enter BEDLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter BEDLOAD Rouse number (W):'};
        name_promptV3='Vertical 3';
        numlines=1;
        onevert=inputdlg(promptV3,name_promptV3,numlines);
        distV3=str2double(onevert(1)); %V3 Distance
        Cz0V3=str2double(onevert(2)); %Base concentration  WL
        W_V3 = str2double(onevert(3)); %Rouse exposant WL
        Cz02V3=str2double(onevert(4)); %Base concentration BL
        W2_V3 = str2double(onevert(5)); %Rouse exposant BL
        
        indV1 = find(abs((round(DMG(:))- distV1))<= 1.5,1, 'first' );
        indV2 = find(abs((round(DMG(:))- distV2))<= 1.5,1, 'first' );
        indV3 = find(abs((round(DMG(:))- distV3))<= 1.5,1, 'first' );
        
        for i = 1:floor((indV1+indV2)/2) %Loop on verticals
            h = max(Depth_Rouse(:,i))+0.5;
            z0 = h-Depth_Rouse(:,i);
            a = 0.05*h;
            for j = 1:numel(Depth_Rouse(:,1)) %Loop on bins
                Conc_fine(j,i) =  Cz0V1.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W_V1;
                Conc_coarse(j,i) =  Cz02V1.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W2_V1;
            end
        end
        
        for i = floor((indV1+indV2)/2)+1:floor((indV3+indV2)/2) %Loop on verticals
            h = max(Depth_Rouse(:,i))+0.5;
            z0 = h-Depth_Rouse(:,i);
            a = 0.05*h;
            for j = 1:numel(Depth_Rouse(:,1)) %Loop on bins
                Conc_fine(j,i) =  Cz0V2.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W_V2;
                Conc_coarse(j,i) =  Cz02V2.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W2_V2;
            end
        end
        
        for i = floor((indV3+indV2)/2)+1:noe %Loop on verticals
            h = max(Depth_Rouse(:,i))+0.5;
            z0 = h-Depth_Rouse(:,i);
            a = 0.05*h;
            for j = 1:numel(Depth_Rouse(:,1)) %Loop on bins
                Conc_fine(j,i) =  Cz0V3.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W_V3;
                Conc_coarse(j,i) =  Cz02V3.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W2_V3;
            end
        end
    case 4
        %% CASE WITH FOUR VERTICALS.
        promptV1={sprintf('Distance from %s', shore),'Enter WASHLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter WASHLOAD Rouse number (W):', 'Enter BEDLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter BEDLOAD Rouse number (W):'};
        name_promptV1='Vertical 1';
        numlines=1;
        onevert=inputdlg(promptV1,name_promptV1,numlines);
        distV1=str2double(onevert(1)); %V1 Distance
        Cz0V1=str2double(onevert(2)); %Base concentration  WL
        W_V1 = str2double(onevert(3)); %Rouse exposant WL
        Cz02V1=str2double(onevert(4)); %Base concentration BL
        W2_V1 = str2double(onevert(5)); %Rouse exposant BL
        %2nd vertical
        promptV2={sprintf('Distance from %s', shore),'Enter WASHLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter WASHLOAD Rouse number (W):', 'Enter BEDLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter BEDLOAD Rouse number (W):'};
        name_promptV2='Vertical 2';
        numlines=1;
        onevert=inputdlg(promptV2,name_promptV2,numlines);
        distV2=str2double(onevert(1)); %V1 Distance
        Cz0V2=str2double(onevert(2)); %Base concentration  WL
        W_V2 = str2double(onevert(3)); %Rouse exposant WL
        Cz02V2=str2double(onevert(4)); %Base concentration BL
        W2_V2 = str2double(onevert(5)); %Rouse exposant BL
        % 3rd vertical
        promptV3={sprintf('Distance from %s', shore),'Enter WASHLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter WASHLOAD Rouse number (W):', 'Enter BEDLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter BEDLOAD Rouse number (W):'};
        name_promptV3='Vertical 3';
        numlines=1;
        onevert=inputdlg(promptV3,name_promptV3,numlines);
        distV3=str2double(onevert(1)); %V3 Distance
        Cz0V3=str2double(onevert(2)); %Base concentration  WL
        W_V3 = str2double(onevert(3)); %Rouse exposant WL
        Cz02V3=str2double(onevert(4)); %Base concentration BL
        W2_V3 = str2double(onevert(5)); %Rouse exposant BL
        %4th vertical
        promptV4={sprintf('Distance from %s', shore),'Enter WASHLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter WASHLOAD Rouse number (W):', 'Enter BEDLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter BEDLOAD Rouse number (W):'};
        name_promptV4='Vertical 4';
        numlines=1;
        onevert=inputdlg(promptV4,name_promptV4,numlines);
        distV4=str2double(onevert(1)); %V4 Distance
        Cz0V4=str2double(onevert(2)); %Base concentration  WL
        W_V4 = str2double(onevert(3)); %Rouse exposant WL
        Cz02V4=str2double(onevert(4)); %Base concentration BL
        W2_V4 = str2double(onevert(5)); %Rouse exposant BL
        
        indV1 = find(abs((round(DMG(:))- distV1))<= 1.5,1, 'first' );
        indV2 = find(abs((round(DMG(:))- distV2))<= 1.5,1, 'first' );
        indV3 = find(abs((round(DMG(:))- distV3))<= 1.5,1, 'first' );
        indV4 = find(abs((round(DMG(:))- distV4))<= 1.5,1, 'first' );
        
        
        for i = 1:floor((indV1+indV2)/2) %Loop on verticals
            h = max(Depth_Rouse(:,i))+0.5;
            z0 = h-Depth_Rouse(:,i);
            a = 0.05*h;
            for j = 1:numel(Depth_Rouse(:,1)) %Loop on bins
                Conc_fine(j,i) =  Cz0V1.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W_V1;
                Conc_coarse(j,i) =  Cz02V1.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W2_V1;
            end
        end
        
        for i = floor((indV1+indV2)/2)+1:floor((indV3+indV2)/2) %Loop on verticals
            h = max(Depth_Rouse(:,i))+0.5;
            z0 = h-Depth_Rouse(:,i);
            a = 0.05*h;
            for j = 1:numel(Depth_Rouse(:,1)) %Loop on bins
                Conc_fine(j,i) =  Cz0V2.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W_V2;
                Conc_coarse(j,i) =  Cz02V2.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W2_V2;
            end
        end
        
        for i = floor((indV3+indV2)/2)+1:floor((indV3+indV4)/2) %Loop on verticals
            h = max(Depth_Rouse(:,i))+0.5;
            z0 = h-Depth_Rouse(:,i);
            a = 0.05*h;
            for j = 1:numel(Depth_Rouse(:,1)) %Loop on bins
                Conc_fine(j,i) =  Cz0V3.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W_V3;
                Conc_coarse(j,i) =  Cz02V3.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W2_V3;
            end
        end
        
        for i = floor((indV3+indV4)/2)+1:noe %Loop on verticals
            h = max(Depth_Rouse(:,i))+0.5;
            z0 = h-Depth_Rouse(:,i);
            a = 0.05*h;
            for j = 1:numel(Depth_Rouse(:,1)) %Loop on bins
                Conc_fine(j,i) =  Cz0V4.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W_V4;
                Conc_coarse(j,i) =  Cz02V4.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W2_V4;
            end
        end
        
        
    case 5
        %% CASE WITH FOUR VERTICALS.
        promptV1={sprintf('Distance from %s', shore),'Enter WASHLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter WASHLOAD Rouse number (W):', 'Enter BEDLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter BEDLOAD Rouse number (W):'};
        name_promptV1='Vertical 1';
        numlines=1;
        onevert=inputdlg(promptV1,name_promptV1,numlines);
        distV1=str2double(onevert(1)); %V1 Distance
        Cz0V1=str2double(onevert(2)); %Base concentration  WL
        W_V1 = str2double(onevert(3)); %Rouse exposant WL
        Cz02V1=str2double(onevert(4)); %Base concentration BL
        W2_V1 = str2double(onevert(5)); %Rouse exposant BL
        %2nd vertical
        promptV2={sprintf('Distance from %s', shore),'Enter WASHLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter WASHLOAD Rouse number (W):', 'Enter BEDLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter BEDLOAD Rouse number (W):'};
        name_promptV2='Vertical 2';
        numlines=1;
        onevert=inputdlg(promptV2,name_promptV2,numlines);
        distV2=str2double(onevert(1)); %V1 Distance
        Cz0V2=str2double(onevert(2)); %Base concentration  WL
        W_V2 = str2double(onevert(3)); %Rouse exposant WL
        Cz02V2=str2double(onevert(4)); %Base concentration BL
        W2_V2 = str2double(onevert(5)); %Rouse exposant BL
        % 3rd vertical
        promptV3={sprintf('Distance from %s', shore),'Enter WASHLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter WASHLOAD Rouse number (W):', 'Enter BEDLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter BEDLOAD Rouse number (W):'};
        name_promptV3='Vertical 3';
        numlines=1;
        onevert=inputdlg(promptV3,name_promptV3,numlines);
        distV3=str2double(onevert(1)); %V3 Distance
        Cz0V3=str2double(onevert(2)); %Base concentration  WL
        W_V3 = str2double(onevert(3)); %Rouse exposant WL
        Cz02V3=str2double(onevert(4)); %Base concentration BL
        W2_V3 = str2double(onevert(5)); %Rouse exposant BL
        %4th vertical
        promptV4={sprintf('Distance from %s', shore),'Enter WASHLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter WASHLOAD Rouse number (W):', 'Enter BEDLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter BEDLOAD Rouse number (W):'};
        name_promptV4='Vertical 4';
        numlines=1;
        onevert=inputdlg(promptV4,name_promptV4,numlines);
        distV4=str2double(onevert(1)); %V4 Distance
        Cz0V4=str2double(onevert(2)); %Base concentration  WL
        W_V4 = str2double(onevert(3)); %Rouse exposant WL
        Cz02V4=str2double(onevert(4)); %Base concentration BL
        W2_V4 = str2double(onevert(5)); %Rouse exposant BL
        %5th vertical
        promptV5={sprintf('Distance from %s', shore),'Enter WASHLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter WASHLOAD Rouse number (W):', 'Enter BEDLOAD bed concentration (Cz0) [mg/L]:',...
            'Enter BEDLOAD Rouse number (W):'};
        name_promptV5='Vertical 5';
        numlines=1;
        onevert=inputdlg(promptV5,name_promptV5,numlines);
        distV5=str2double(onevert(1)); %V5 Distance
        Cz0V5=str2double(onevert(2)); %Base concentration  WL
        W_V5 = str2double(onevert(3)); %Rouse exposant WL
        Cz02V5=str2double(onevert(4)); %Base concentration BL
        W2_V5 = str2double(onevert(5)); %Rouse exposant BL
        
        indV1 = find(abs((round(DMG(:))- distV1))<= 1.5,1, 'first' );
        indV2 = find(abs((round(DMG(:))- distV2))<= 1.5,1, 'first' );
        indV3 = find(abs((round(DMG(:))- distV3))<= 1.5,1, 'first' );
        indV4 = find(abs((round(DMG(:))- distV4))<= 1.5,1, 'first' );
        indV5 = find(abs((round(DMG(:))- distV5))<= 1.5,1, 'first' );
        for i = 1:floor((indV1+indV2)/2) %Loop on verticals
            h = max(Depth_Rouse(:,i))+0.5;
            z0 = h-Depth_Rouse(:,i);
            a = 0.05*h;
            for j = 1:numel(Depth_Rouse(:,1)) %Loop on bins
                Conc_fine(j,i) =  Cz0V1.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W_V1;
                Conc_coarse(j,i) =  Cz02V1.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W2_V1;
            end
        end
        
        for i = floor((indV1+indV2)/2)+1:floor((indV3+indV2)/2) %Loop on verticals
            h = max(Depth_Rouse(:,i))+0.5;
            z0 = h-Depth_Rouse(:,i);
            a = 0.05*h;
            for j = 1:numel(Depth_Rouse(:,1)) %Loop on bins
                Conc_fine(j,i) =  Cz0V2.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W_V2;
                Conc_coarse(j,i) =  Cz02V2.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W2_V2;
            end
        end
        
        for i = floor((indV3+indV2)/2)+1:floor((indV3+indV4)/2) %Loop on verticals
            h = max(Depth_Rouse(:,i))+0.5;
            z0 = h-Depth_Rouse(:,i);
            a = 0.05*h;
            for j = 1:numel(Depth_Rouse(:,1)) %Loop on bins
                Conc_fine(j,i) =  Cz0V3.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W_V3;
                Conc_coarse(j,i) =  Cz02V3.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W2_V3;
            end
        end
        
        for i = floor((indV3+indV4)/2)+1:floor((indV5+indV4)/2) %Loop on verticals
            h = max(Depth_Rouse(:,i))+0.5;
            z0 = h-Depth_Rouse(:,i);
            a = 0.05*h;
            for j = 1:numel(Depth_Rouse(:,1)) %Loop on bins
                Conc_fine(j,i) =  Cz0V4.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W_V4;
                Conc_coarse(j,i) =  Cz02V4.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W2_V4;
            end
        end
        
        for i = floor((indV5+indV4)/2)+1:noe %Loop on verticals
            h = max(Depth_Rouse(:,i))+0.5;
            z0 = h-Depth_Rouse(:,i);
            a = 0.05*h;
            for j = 1:numel(Depth_Rouse(:,1)) %Loop on bins
                Conc_fine(j,i) =  Cz0V5.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W_V5;
                Conc_coarse(j,i) =  Cz02V5.*(((h-z0(j))./z0(j)).*(a./(h-a))).^W2_V5;
            end
        end
end
%% Full concentrations
Conc_full = Conc_coarse+Conc_fine;

for i = 1:noe
    Conc_full(nlin_depth(i):end,i)=NaN; % Convert back unmeasureable values to Nan (below the bed)
    Conc_fine(nlin_depth(i):end,i)=NaN; % Convert back unmeasureable values to Nan (below the bed)
    Conc_coarse(nlin_depth(i):end,i)=NaN; % Convert back unmeasureable values to Nan (below the bed)
end
%% Velocity
vel2recons_log = log10(vel2recons);
velfull_log = inpaint_nans(vel2recons_log,4);% Perform on log of velocities as the folow a logarthmic profile on the verticals
velfull = 10.^(velfull_log);
for i = 1:noe
    velfull(nlin_depth(i):end,i)=NaN; % Convert back unmeasureable values to Nan (below the bed)
end

%% Q
Q2recons_log = log10(Q2recons);
%Q2recons_log = Q2recons;
Qfull_log = nan_interp(Q2recons_log,2);% Perform on log of Qocities as the folow a logarthmic profile on the verticals
Qfull = 10.^(Qfull_log);
%Qfull = Qfull_log;
for i = 1:noe
    Qfull(nlin_depth(i):end,i)=NaN;
end

%% FLUX CALCULATION


full_flux = Qfull.*(Conc_full./1000); %kg/s
fine_flux = Qfull.*(Conc_fine./1000); %kg/s
coarse_flux = Qfull.*(Conc_coarse./1000); %kg/s
%fine_flux = Q2recons.*Conc_fine;



% VELOCITY
figure;
pcolor(DMG,Depth_binsTot(:,1),velfull)
hold on
shading interp
colormap(jet)
t=colorbar;
xlabel('Distance [m]');
set(gca,'YDir','reverse');
set(get(t,'ylabel'),'String', 'Velocity [m/s]');
ylabel('Depth [m]');
title(sprintf('Extrapolated velocities at %s',...
    namedate),'FontSize', 10, 'FontWeight','bold')
hold on;
plot(DMG, DepthBave, 'color', 'black', 'LineWidth',1.5)
hold on;
plot(DMG, endDepth_Plot, 'color', 'black', 'LineWidth',0.75)
hold off;

% FLOW
figure;
pcolor(DMG,Depth_binsTot(:,1),Qfull)
hold on
shading interp
colormap(jet)
t=colorbar;
xlabel('Distance [m]');
set(gca,'YDir','reverse');
set(get(t,'ylabel'),'String', 'Flow [m3/s]');
ylabel('Depth [m]');
title(sprintf('Extrapolated discharge at %s',...
    namedate),'FontSize', 10, 'FontWeight','bold')
hold on;
plot(DMG, DepthBave, 'color', 'black', 'LineWidth',1.5)
hold on;
plot(DMG, endDepth_Plot, 'color', 'black', 'LineWidth',0.75)
hold off;

% CONCENTRATION
%full
figure;
subplot(3,1,1)
pcolor(DMG,Depth_binsTot(:,1),Conc_full)
hold on
shading interp
colormap(jet)
t=colorbar;
xlabel('Distance [m]');
set(gca,'YDir','reverse');
set(get(t,'ylabel'),'String', 'Sediment concentration [mg/L]');
ylabel('Depth [m]');
title(sprintf('Extrapolated sediment concentration at %s',...
    namedate),'FontSize', 10, 'FontWeight','bold')
hold on;
plot(DMG, DepthBave, 'color', 'black', 'LineWidth',1.5)
hold on;
plot(DMG, endDepth_Plot, 'color', 'black', 'LineWidth',0.75)
hold off;
%fine
subplot(3,1,2)
pcolor(DMG,Depth_binsTot(:,1),Conc_fine)
hold on
shading interp
colormap(jet)
t=colorbar;
xlabel('Distance [m]');
set(gca,'YDir','reverse');
set(get(t,'ylabel'),'String', 'Sediment concentration [mg/L]');
ylabel('Depth [m]');
title(sprintf('Extrapolated fine sediment concentration at %s',...
    namedate),'FontSize', 10, 'FontWeight','bold')
hold on;
plot(DMG, DepthBave, 'color', 'black', 'LineWidth',1.5)
hold on;
plot(DMG, endDepth_Plot, 'color', 'black', 'LineWidth',0.75)
hold off;
%Coarse
subplot(3,1,3)
pcolor(DMG,Depth_binsTot(:,1),Conc_coarse)
hold on
shading interp
colormap(jet)
t=colorbar;
xlabel('Distance [m]');
set(gca,'YDir','reverse');
set(get(t,'ylabel'),'String', 'Sediment concentration [mg/L]');
ylabel('Depth [m]');
title(sprintf('Extrapolated coarse sediment concentration at %s',...
    namedate),'FontSize', 10, 'FontWeight','bold')
hold on;
plot(DMG, DepthBave, 'color', 'black', 'LineWidth',1.5)
hold on;
plot(DMG, endDepth_Plot, 'color', 'black', 'LineWidth',0.75)
hold off;

% FLUX
%full
figure;
subplot(3,1,1)
pcolor(DMG,Depth_binsTot(:,1),full_flux)
hold on
shading interp
colormap(jet(512))
t=colorbar;
xlabel('Distance [m]');
set(gca,'YDir','reverse');
set(get(t,'ylabel'),'String', 'Sediment flux [kg/s]');
ylabel('Depth [m]');
title(sprintf('Calculated total sediment flux at %s',...
    namedate),'FontSize', 10, 'FontWeight','bold')
hold on;
plot(DMG, DepthBave, 'color', 'black', 'LineWidth',1.5)
hold on;
plot(DMG, endDepth_Plot, 'color', 'black', 'LineWidth',0.75)
hold off;
%fine
subplot(3,1,2)
pcolor(DMG,Depth_binsTot(:,1),fine_flux)
hold on
shading interp
colormap(jet(512))
t=colorbar;
xlabel('Distance [m]');
set(gca,'YDir','reverse');
set(get(t,'ylabel'),'String', 'Sediment flux [kg/s]');
ylabel('Depth [m]');
title(sprintf('Calculated fine sediment flux at %s',...
    namedate),'FontSize', 10, 'FontWeight','bold')
hold on;
plot(DMG, DepthBave, 'color', 'black', 'LineWidth',1.5)
hold on;
plot(DMG, endDepth_Plot, 'color', 'black', 'LineWidth',0.75)
hold off;
%coarse
subplot(3,1,3)
pcolor(DMG,Depth_binsTot(:,1),coarse_flux)
hold on
shading interp
colormap(jet(512))
t=colorbar;
xlabel('Distance [m]');
set(gca,'YDir','reverse');
set(get(t,'ylabel'),'String', 'Sediment flux [kg/s]');
ylabel('Depth [m]');
title(sprintf('Calculated coarse sediment flux at %s',...
    namedate),'FontSize', 10, 'FontWeight','bold')
hold on;
plot(DMG, DepthBave, 'color', 'black', 'LineWidth',1.5)
hold on;
plot(DMG, endDepth_Plot, 'color', 'black', 'LineWidth',0.75)
hold off;
%% Diff Q tot RDI and interp
Q_tot= A.Q.top+A.Q.meas+A.Q.start+A.Q.end;
diffRDI_interp = (nansum(nansum(Qfull))-max(Q_tot))/nansum(nansum(Qfull))*100;
totQinterp = nansum(nansum(Qfull));
fprintf('\nTotal flow: %6.3f [m3/s]\n',totQinterp);
fprintf('\nDifference of flow value between the RDI derived value and the interpolated value: %6.1f percent \n',diffRDI_interp);


%% Disp tot flux fine
totflux = double(nansum(nansum(full_flux)));
totfluxfine = double(nansum(nansum(fine_flux)));
totfluxcoarse = double(nansum(nansum(coarse_flux)));
totflux_day = (totflux/1000)*24*60*60;
totflux_month = totflux_day*30;
meanconc_fine = double(nanmean(nanmean(Conc_fine)));
meanconc_coarse = double(nanmean(nanmean(Conc_coarse)));
meanconc_full = double(nanmean(nanmean(Conc_full)));
fprintf('\nMean concentration of fines: %6.2f [mg/L]\n', meanconc_fine);
fprintf('\nMean concentration of coarses: %6.2f [mg/L]\n', meanconc_coarse);
fprintf('\nMean total concentration: %6.2f [mg/L]\n', meanconc_full);
fprintf('\nFine flux: %6.2f [kg/s] \n',totfluxfine);
fprintf('\nCoarse flux: %6.2f [kg/s] \n',totfluxcoarse);
fprintf('\nTotal flux: %6.2f [kg/s] \n',totflux);
fprintf('\nTotal flux: %6.2f [t/day] \n',totflux_day);
fprintf('\nTotal flux: %6.2f [t/month] \n',totflux_month);

%% Part to see if the user wishes to compare the results with measured concentrations of a vertical
correcSelec = questdlg('Do you want to compare your results with measured concentrations on a vertical?', ...
    'Comparison with in-situ vertical data', ...
    'Yes','No','Yes');
switch correcSelec
    case 'Yes'
        
        
        
        %Select concentration file to compare the results:
        [Conc_file1, pathname_conc] = ...
            uigetfile({'*.txt'},'Select in-situ concentrations (N.B format = [depth Conc])');
        Conc_file2 = [pathname_conc Conc_file1];
        Conc_file = load(Conc_file2);
        Conc_depth = Conc_file(:,1);
        TotC=Conc_file(:,2); % Total concentration [mg/l]
        FineC=Conc_file(:,3); % Fines sediment concentration [mg/l]
        CoarseC=Conc_file(:,4); % Coarse sediment concentration [mg/l]
        %PercentFines = Conc_file(:,5); % % Fines
        %PercentCoarse = Conc_file(:,6); % % Coarse
        % Ask user to insert distance from ban and number of verticals to analyse
        promptCheck={'Distance between measurement and  shore','Number of ensembles/verticals to average:'};
        name_promptCheck='Comparison of results';
        numlines=1;
        Checkans= {'1200','100'};
        options.Resize='on';
        options.WindowStyle='normal';
        options.Interpreter='tex';
        Checkparam=inputdlg(promptCheck,name_promptCheck,numlines,Checkans, options);
        dist=str2double(Checkparam(1)); % median diameter (m)
        nb_vert=str2double(Checkparam(2)); % standard deviation
        
        %find the location of the closest vertical corresponding to the distance
        indVert = find(abs((round(DMG(:))- dist))<= 1.5,1, 'first' );
        
        %Get the indices of the verticals to average
        indAV = indVert-ceil(nb_vert/2)+1:1:indVert+ceil(nb_vert/2);
        
        %isolate number of verticals for the velocity magnitude to average
        vel2av = velfull(:,indAV);
        meanV = nanmean(vel2av,2);
        depth2av = DepthBave(indAV);
        avDepth = nanmean(depth2av);
        % IntB1_2av = Intensity(:,indAV);
        % avIntB1 = nanmean(IntB1_2av,2);
        DepthB12av=DepthBave(indAV);
        avDepthB1 = nanmean(DepthB12av);
        Rangesvec2av = A.Wat.binDepth(:,indAV);
        avRangesvec = nanmean(Rangesvec2av,2);
        Conc_fine2av = Conc_fine(:,indAV);
        Conc_coarse2av = Conc_coarse(:,indAV);
        avConc_coarse = nanmean(Conc_coarse2av,2);
        avConc_fine = nanmean(Conc_fine2av,2);
        Depthfull2av = Depth_binsTot(:,indAV);
        
        
        %% Plots of comparison of estimated and measured concentrations
        
        % Plot velocity as a function of depth
        figure;
        for i = 1:numel(indAV)
            plot(vel2av(:,i),-Depthfull2av(:,i))
            hold on
            plot(meanV,-Depthfull2av(:,ceil(numel(indAV)/2)), 'LineWidth',2, 'Color','r')
        end
        
        ylabel('Depth [m]')
        xlabel('Velocity [m/s]')
        legend('Ensemble velocities', 'Mean velocity')
        title('Velocities of the ensembles and mean velocity for the considered vertical')
        
        % Plot velocity as a function of depth normalized by mean depth
        figure;
        
        for i = 1:numel(indAV)
            plot(vel2av(:,i),1-(Depthfull2av(:,i)./avDepthB1))
            hold on
            plot(meanV,1-(Depthfull2av(:,ceil(numel(indAV)/2))./avDepthB1), 'LineWidth',2, 'Color','r')
            
        end
        ylim('manual')
        ylim = ([0 1]);
        ylabel('Dimensionless elevation above the bed')
        xlabel('Velocity [m/s]')
        legend('Ensemble velocities', 'Mean velocity')
        %title('Velocities of the ensembles and mean velocity for the considered vertical')
        title('Velocities of the ensembles and mean velocity for the considered vertical ',...
            'FontSize', 11, 'FontWeight','bold')
        
        
        % Plot fine sediment concentration as a function of depth normalized by mean depth
        figure;
        
        plot(avConc_fine,1-(Depthfull2av(:,ceil(numel(indAV)/2))./avDepthB1), 'LineWidth',2, 'Color','r')
        hold on
        plot(FineC, 1-(Conc_depth./avDepthB1), 'LineStyle','none','Marker', 'o','Color', 'black', 'MarkerSize',3, 'MarkerFaceColor','black')
        %ylim('manual')
        %ylim = ([0 1]);
        ylabel('Dimensionless elevation above the bed')
        xlabel('Concentration [kg/m3]')
        legend('Modelled washload', 'Measured washload')
        %title('Velocities of the ensembles and mean velocity for the considered vertical')
        title('Measured and observed washload concentrations for the considered vertical',...
            'FontSize',9, 'FontWeight','bold')
        figure;
        
        plot(avConc_coarse,1-(Depthfull2av(:,ceil(numel(indAV)/2))./avDepthB1), 'LineWidth',2, 'Color','r')
        hold on
        plot(CoarseC, 1-(Conc_depth./avDepthB1), 'LineStyle','none','Marker', 'o','Color', 'black', 'MarkerSize',3, 'MarkerFaceColor','black')
        %ylim('manual')
        %ylim = ([0 1]);
        ylabel('Dimensionless elevation above the bed')
        xlabel('Concentration [kg/m3]')
        legend('Modelled bedload', 'Measured bedload')
        %title('Velocities of the ensembles and mean velocity for the considered vertical')
        title('Measured and observed bedload concentrations for the considered vertical',...
            'FontSize',9, 'FontWeight','bold')
        
        figure;
        
        plot((avConc_fine+avConc_coarse),1-(Depthfull2av(:,ceil(numel(indAV)/2))./avDepthB1), 'LineWidth',2, 'Color','r')
        hold on
        plot((CoarseC+FineC), 1-(Conc_depth./avDepthB1), 'LineStyle','none','Marker', 'o','Color', 'black', 'MarkerSize',3, 'MarkerFaceColor','black')
        %ylim('manual')
        %ylim = ([0 1]);
        ylabel('Dimensionless elevation above the bed')
        xlabel('Concentration [kg/m3]')
        legend('Modelled total load', 'Measured total load')
        %title('Velocities of the ensembles and mean velocity for the considered vertical')
        title('Measured and observed total load concentrations for the considered vertical',...
            'FontSize',9, 'FontWeight','bold')
        
    case 'No'
        disp('End of processing')
end
