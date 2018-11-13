function [A]=readADCPfile(fullName,screenData,ignoreBS);
% tfile reads the data from an RDI ASCII output file and puts the
% data in a Matlab data structure with major groups of:
% Sup - supporing data
% Wat - water data
% Nav - navigation data including GPS.
% Sensor - Sensor data
% Q - discharge related data
%
% The data can be screened (screenData=1) so that invalid data are set to
% nan or data reflect strictly the ASCII output file (screenData=0). If
% screenData=0 then the data reflect the ASCII file and -32768 and other
% RDI defined values are left in the data structure. If screenData=1 the
% RDI defined values are trapped and set to nan for many variables.
%
% WinRiver II will sometimes put a $ instead of the correct values for
% intensity of backscatter. Setting ignoreBS=1 will skip decoding the
% intensity of backscatter data to compensate for this bug in WR2.
%
% David S. Mueller, USGS, Office of Surface Water
% Frank L. Engel, USGS, Illinois Water Science Center

%% Initial Scan of File
%  Initial scan required to preallocate arrays. All ensembles must be
%  scanned because RiverRay has variable number of bins.

% Turn off LaTex interpreter to avoid subscripts in filenames containing _
set(0, 'DefaulttextInterpreter', 'none');
fid=fopen(fullName);
% clc
% disp(['Scanning Data File: ' fullName]);
idx=find(fullName=='\',1,'last');
fileName=fullName(idx:end);

% Display waitbar
waitmessage=['Reading ' fileName];
hwait=waitbar(0,waitmessage);

% Scan Fixed Leader.
lineIn=fgetl(fid);
lineIn=fgetl(fid);
lineIn=fgetl(fid);

% Count lines in file.
k=1;
fileEnd=0;

% Loop required to determine number of variable leader lines due to bug
% in ASCII output from WinRiver II for reference set to NONE.
check=0;
linecount=0;
while check==0
    lineIn=fgetl(fid);
    idxBT=strfind(lineIn,'BT');
    idxGGA=strfind(lineIn,'GGA');
    idxVTG=strfind(lineIn,'VTG');
    idxNone=strfind(lineIn,'NONE');
    check=nansum([idxBT idxGGA idxVTG idxNone]);
    linecount=linecount+1;
end

% Read number of bins in 1st ensemble
[bins(k),~,~,~,~,~]=strread(lineIn,'%f %s %s %s %f %f',1);

% Skip bin data
dummy=textscan(fid, '%s %*[^\n]',bins(k));

% Set number of leader lines to get to number of bins in variable
% leader
leaderlines=linecount-1;

% Begin loop to determin number of ensembles and maximum number of bins
% for preallocating arrays. Looping through the entire file is required
% because RiverRay data have variable numbers of bins in each ensemble.
while fileEnd==0
    dummy=textscan(fid, '%s %*[^\n]',leaderlines);
    if length(dummy{1})>1
        k=k+1;
        bins(k)=cell2mat(textscan(fid, '%f %*[^\n]',1));
        dummy=textscan(fid, '%s %*[^\n]',bins(k));
    end
    fileEnd=feof(fid);
end

% Update waitbar
waitbar(0.1);
%% Complete scan of input file.
% Close input file and report that the initial scan is completed.
fclose(fid);
% disp('Scan Complete');

%% Initialize Data Structure
% Preallocates arrays based on information from initial scan.

% This is the number of ensembles actually contained in the input file.
noe=k;
bins=max(bins);

% Initialize Data Structure.
Sup=struct( 'absorption_dbpm',nan(noe,1),...
    'bins',nan(noe,1),...
    'binSize_cm',nan(1),...
    'nBins',nan(1),...
    'blank_cm',nan(1),...
    'draft_cm',nan(1),...
    'ensNo',nan(noe,1),...
    'nPings',nan(1),...
    'noEnsInSeg',nan(noe,1),...
    'noe',nan(1),...
    'note1',blanks(80),...
    'note2',blanks(80),...
    'intScaleFact_dbpcnt',nan(noe,1),...
    'intUnits',repmat(blanks(5),noe,1),...
    'vRef',repmat(blanks(4),noe,1),...
    'wm',nan(1),...
    'units',repmat(blanks(2),noe,1),...
    'year',nan(noe,1),...
    'month',nan(noe,1),...
    'day',nan(noe,1),...
    'hour',nan(noe,1),...
    'minute',nan(noe,1),...
    'second',nan(noe,1),...
    'sec100',nan(noe,1),...
    'timeElapsed_sec',nan(noe,1),...
    'timeDelta_sec100',nan(1));

Wat=struct( 'binDepth',nan(bins,noe),...
    'backscatter',nan(bins,noe,4),...
    'vDir',nan(bins,noe),...
    'vMag',nan(bins,noe),...
    'vEast',nan(bins,noe),...
    'vError',nan(bins,noe),...
    'vNorth',nan(bins,noe),...
    'vVert',nan(bins,noe),...
    'percentGood',nan(bins,noe));

Nav=struct( 'bvEast',nan(noe,1),...
    'bvError',nan(noe,1),...
    'bvNorth',nan(noe,1),...
    'bvVert',nan(noe,1),...
    'depth',nan(noe,4),...
    'dsDepth',nan(noe,1),...
    'dmg',nan(noe,1),...
    'length',nan(noe,1),...
    'totDistEast',nan(noe,1),...
    'totDistNorth',nan(noe,1),...
    'altitude',nan(noe,1),...
    'altitudeChng',nan(noe,1),...
    'gpsTotDist',nan(noe,1),...
    'gpsVariable',nan(noe,1),...
    'gpsVeast',nan(noe,1),...
    'gpsVnorth',nan(noe,1),...
    'lat_deg',nan(noe,1),...
    'long_deg',nan(noe,1),...
    'nSats',nan(noe,1),...
    'hdop',nan(noe,1));

Sensor=struct(  'pitch_deg',nan(noe,1),...
    'roll_deg',nan(noe,1),...
    'heading_deg',nan(noe,1),...
    'temp_degC',nan(noe,1));

Q=struct(   'endDepth',nan(noe,1),...
    'endDist',nan(noe,1),...
    'bot',nan(noe,1),...
    'end',nan(noe,1),...
    'meas',nan(noe,1),...
    'start',nan(noe,1),...
    'top',nan(noe,1),...
    'unit',nan(bins,noe),...
    'startDepth',nan(noe,1),...
    'startDist',nan(noe,1));

Sup.noe=noe;

%% Read File and Store Data
% All data are read and stored in preallocated data structures

% Reopen File for Reading
fid=fopen(fullName);
% disp('Reading Data File');

% Read Fixed Leader
Sup.note1=fgetl(fid);
Sup.note2=fgetl(fid);
C = textscan(fid,'%u %u %u %u %u %u %u',1);
[...
    Sup.binSize_cm,...
    Sup.blank_cm,...
    Sup.draft_cm,...
    Sup.nBins,...
    Sup.nPings,...
    Sup.timeDelta_sec100,...
    Sup.wm] = deal(C{:});

% Read Variable Leader
waitstep = floor(noe/100);
for n = 1:noe;
    % Update the waitbar only on whole percents
    if ~mod(n, waitstep) || n==noe
        waitbar(0.1+n/noe);
    end
    C = textscan(fid,'%u %u %u %u %u %f %f %u %u %f %f %f %f',1);
    [...
        Sup.year(n),...
        Sup.month(n),...
        Sup.day(n),...
        Sup.hour(n),...
        Sup.minute(n),...
        Sup.second(n),...
        Sup.sec100(n),...
        Sup.ensNo(n),...
        Sup.noEnsInSeg(n),...
        Sensor.pitch_deg(n),...
        Sensor.roll_deg(n),...
        Sensor.heading_deg(n),...
        Sensor.temp_degC(n)] = deal(C{:});
    
    % Required logic to account for bug in ASCII output when the
    % reference is set to NONE.
    C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f',1);
    if sum(isnan([C{:}]))==0
        [...
            Nav.bvEast(n),...
            Nav.bvNorth(n),...
            Nav.bvVert(n),...
            Nav.bvError(n),...
            Nav.dsDepth(n),...
            Nav.altitude(n),...
            Nav.altitudeChng(n),...
            Nav.gpsVariable(n),...
            Nav.depth(n,1),...
            Nav.depth(n,2),...
            Nav.depth(n,3),...
            Nav.depth(n,4)] = deal(C{:});
        
        C = textscan(fid,'%f %f %f %f %f',1);
        [...
            Nav.length(n),...
            Sup.timeElapsed_sec(n),...
            Nav.totDistNorth(n),...
            Nav.totDistEast(n),...
            Nav.dmg(n)] = deal(C{:});
        
        C = textscan(fid,'%f %f %f %f %f',1);
        [...
            Nav.lat_deg(n),...
            Nav.long_deg(n),...
            Nav.gpsVeast(n),...
            Nav.gpsVnorth(n),...
            Nav.gpsTotDist(n)] = deal(C{:});
    else
        [...
            Nav.dsDepth(n),...
            Nav.altitude(n),...
            Nav.altitudeChng(n),...
            Nav.gpsVariable(n),...
            Nav.depth(n,1),...
            Nav.depth(n,2),...
            Nav.depth(n,3),...
            Nav.depth(n,4)]=deal(C{1:8});
        
        C = textscan(fid,'%f %f %f %f %f %f',1);
        [...
            dummy,...
            Nav.lat_deg(n),...
            Nav.long_deg(n),...
            Nav.gpsVeast(n),...
            Nav.gpsVnorth(n),...
            Nav.gpsTotDist(n)] = deal(C{:});
    end
    
    % Extract HDOP and number of satellites from gpsVariable
    if Nav.gpsVariable(n)>0
        Nav.hdop(n)=floor(Nav.gpsVariable(n))./10;  %/10 added 3-12-10 by PRJ (according to TRDI WRII manual)
        Nav.nSats(n)=(Nav.gpsVariable(n)-Nav.hdop(n).*10).*100;
    end;
    
    C = textscan(fid,'%f %f %f %f %f %f %f %f %f',1);
    [...
        Q.meas(n),...
        Q.top(n),...
        Q.bot(n),...
        Q.start(n),...
        Q.startDist(n),...
        Q.end(n),...
        Q.endDist(n),...
        Q.startDepth(n),...
        Q.endDepth(n)] = deal(C{:});
    
    C = textscan(fid,'%f %s %s %s %f %f',1);
    [...
        Sup.bins(n),...
        Sup.units,...
        Sup.vRef,...
        Sup.intUnits,...
        Sup.intScaleFact_dbpcnt(n),...
        Sup.absorption_dbpm(n)]= deal(C{:});
    
    % Read Profile Data.
    
    
    % Logic to account for WR2 bug that puts $ in for intensity or
    % backscatter.
    if ignoreBS==0
        C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f',Sup.bins(n));
        [...
            Wat.binDepth(1:Sup.bins(n),n),...
            Wat.vMag(1:Sup.bins(n),n),...
            Wat.vDir(1:Sup.bins(n),n),...
            Wat.vEast(1:Sup.bins(n),n),...
            Wat.vNorth(1:Sup.bins(n),n),...
            Wat.vVert(1:Sup.bins(n),n),...
            Wat.vError(1:Sup.bins(n),n),...
            Wat.backscatter(1:Sup.bins(n),n,1),...
            Wat.backscatter(1:Sup.bins(n),n,2),...
            Wat.backscatter(1:Sup.bins(n),n,3),...
            Wat.backscatter(1:Sup.bins(n),n,4),...
            Wat.percentGood(1:Sup.bins(n),n),...
            Q.unit(1:Sup.bins(n),n)] = deal(C{:});
    else
        C = textscan(fid,'%f %f %f %f %f %f %f %s %s %s %s %f %f',Sup.bins(n));
        [...
            Wat.binDepth(1:Sup.bins(n),n),...
            Wat.vMag(1:Sup.bins(n),n),...
            Wat.vDir(1:Sup.bins(n),n),...
            Wat.vEast(1:Sup.bins(n),n),...
            Wat.vNorth(1:Sup.bins(n),n),...
            Wat.vVert(1:Sup.bins(n),n),...
            Wat.vError(1:Sup.bins(n),n),...
            dummy,...
            dummy,...
            dummy,...
            dummy,...
            Wat.percentGood(1:Sup.bins(n),n),...
            Q.unit(1:Sup.bins(n),n)] = deal(C{:});
    end
    
end % for loop through ensembles

% Close File
fclose(fid);
% disp('Data Input Complete');

%% Screen Data
% Screens invalid data replacing them with nan

% Apply data screening if requested.
if screenData==1
    
    Wat.vNorth(Wat.vNorth==-32768)=nan;
    Wat.vEast(Wat.vEast==-32768)=nan;
    Wat.vError(Wat.vError==-32768)=nan;
    Wat.vVert(Wat.vVert==-32768)=nan;
    Wat.vMag(Wat.vMag==-32768)=nan;
    Wat.vDir(Wat.vDir==-32768)=nan;
    Wat.percentGood(Wat.percentGood==-32768)=nan;
    Wat.backscatter(Wat.backscatter==-32768)=nan;
    
    
    Nav.bvNorth(Nav.bvNorth==-32768)=nan;
    Nav.bvEast(Nav.bvEast==-32768)=nan;
    Nav.bvError(Wat.vError==-32768)=nan;
    Nav.bvVert(Nav.bvVert==-32768)=nan;
    Nav.depth(Nav.depth==0)=nan;
    Nav.lat_deg(Nav.lat_deg==30000)=nan;
    Nav.long_deg(Nav.long_deg==30000)=nan;
    Nav.gpsVnorth(Nav.gpsVnorth==-32768)=nan;
    Nav.gpsVeast(Nav.gpsVeast==-32768)=nan;
    Q.unit(Q.unit==2147483647)=nan;
end;

% Assign Data to One Structure
A.Sup=Sup;
A.Wat=Wat;
A.Nav=Nav;
A.Sensor=Sensor;
A.Q=Q;

% Close waitbar
close(hwait);
