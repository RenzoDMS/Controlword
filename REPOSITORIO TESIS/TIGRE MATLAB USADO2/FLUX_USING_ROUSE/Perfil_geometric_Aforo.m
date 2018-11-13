0%%----------------------------------------------------------------------
% récupération de la matrice des vitesses à partir du fichier initial
% consituté de la matrice des profondeurs + matrice des vitesses agrégées
% l'une après l'autre.
%----------------------------------------------------------------------
clear all;
close all;
clc;
%the folder where the files to be read, the script and the functions are.
[filename, pathname_ascii] = ...
    uigetfile({'*.txt'},'Select classic ADCP ASCII output from WRII');

fullName=[pathname_ascii filename]; %Enter here name of classic ASCII output file

ignoreBS=0;
screenData=1;
A=readADCPfile(fullName,screenData,ignoreBS);
% from WinRiver you wish to read
stationname = inputdlg('Please enter the name of the station:', 'Station');
stationname = num2str(stationname{1});
yr = num2str(A.Sup.year(1));
month = num2str(A.Sup.month(1));
day = num2str(A.Sup.day(1));
date = [day,'-', month,'-', yr];
date = num2str(date);
namedate = [stationname, ' on the ', date];



%% Extraction of ADCP data for processing

velocity_temp = A.Wat.vMag./100; %velocity magnitude

% Total distance made good
distance = A.Nav.dmg;

prof_temp = A.Wat.binDepth;

% initialisation des constantes

num_figure = 0; % gestion des figures


%vel2recons= [MatrixUpperPart;velMag ];

% affichage du champ de vitesses brut
%figure (num_figure);
%imagesc(velocity_temp);

%--------------------------------------------------------------------------
%rajout des cellules de surfaces manquantes
%--------------------------------------------------------------------------
%calcul du nombre de lignes à rajouter en surface
hauteur_cellule = prof_temp(2,1)-prof_temp(1,1);
nb_cel_ajout = floor(prof_temp(1)/hauteur_cellule);
noe=A.Sup.noe; % number of ensembles
[nbligtemp nbcol] = size(prof_temp);
nblig = nb_cel_ajout + nbligtemp;
%création des matrices avec les lignes à rajouter en plus
prof = ones(nblig, nbcol);
velocity = ones(nblig, nbcol);

%calcul des profondeurs sur les cellules supérieures
profajout = ones( nb_cel_ajout,1);
base = prof_temp(1,1) - nb_cel_ajout*hauteur_cellule;
for i= 1:nb_cel_ajout
    profajout(i) = base +(i-1)*hauteur_cellule;
end

% remplissage des matrices
for i= 1:nbcol
    prof((nb_cel_ajout+1):nblig,i) = prof_temp(:,i);
    velocity((nb_cel_ajout+1):nblig,i) = velocity_temp(:,i);
end

for i=1:nbcol
    prof(1:nb_cel_ajout,i) = profajout(:);
    velocity(1:nb_cel_ajout,i) = velocity_temp(1,i);
end

% affichage du champ de vitesses brut
num_figure = num_figure+1;
figure(num_figure);
%imagesc(velocity);
distance_temp=distance;
imagesc(distance_temp,prof_temp(:,1),velocity_temp);
title(namedate,'FontName','Arial','FontSize', 14);
ylabel('depth (m)','FontName', 'Arial','FontSize', 14);
xlabel('distance (m)','FontName', 'Arial','FontSize', 14);
colorbar
% élimination des variables auxiliaires
%clear data;
%clear aux_vel;
%clear pos_max_deep;
%clear prof_temp;
%clear velocity_temp;

%%
%--------------------------------------------------------------------
% création du vecteur de profondeur max où une vitesse a été enregistrée
% pour chaque colonne
%---------------------------------------------------------------------
DepthBave = mean(A.Nav.depth,2);
BinSize = double(A.Sup.binSize_cm)/100;
Depth_bins = A.Wat.binDepth;
indepthBave=find(DepthBave<=0);
DepthBave(indepthBave)=NaN;
DepthBave = inpaint_nans(DepthBave,2);
nlin_depth = nan(size(DepthBave,1),1);
ncol_depth = nan(size(DepthBave,1),1);
for i =1:noe
    [nlin_depth(i), ncol_depth(i)]= find(Depth_bins(:,i)<=DepthBave(i) + BinSize & Depth_bins(:,i)>=DepthBave(i)-BinSize, 1,'Last'); %N.B. A.Q.endDepth is the last measured value
    Vtemp0=[nlin_depth, ncol_depth];
end

velocity_log = log10(velocity);
velfull_log = inpaint_nans(velocity_log,4);% Perform on log of velocities as the folow a logarthmic profile on the verticals

velocity = 10.^(velfull_log);


for i = 1:noe
    velocity(nlin_depth(i):end,i)=NaN; % Convert back unmeasureable values to Nan (below the bed)
end


t = size(velocity);
nbcol = t(:,2); % nombre de verticales à traiter 
nblig = t(:,1); % nombre de cellules max ADCP par ensemble
%idp=[]
pmax=ones(1,nbcol);
d=1;              
% création du vecteur de profondeur max avec valeur de vitesse
velnan=isnan(velocity);
for i=d:nbcol,
     pmax(i)=find(velnan(:,i)==0,1,'last');
    % eval(['pmax(' num2str(i) ')=id;']);
end

% création du vecteur de profondeur min sans valeur de vitesse
pmin=ones(1,nbcol);
for i=d:nbcol,            
   pmin(i)=find(velnan(:,i)==1,1,'first');
   % eval(['pmin(' num2str(i) ')=id;']);
end

%% Get the maximum depth in every vertical
for i=1:nbcol
    Yfondo(1,i) =prof(pmax(i),i) ;   
end
Yfondo=Yfondo.';
Yfondo=0-Yfondo;
distance=[0;distance;distance(nbcol)];
Yfondo=[0;Yfondo;0];
XY=[distance Yfondo];
xlswrite('CoordenadasSección.xlsx', XY , 'Ymax-Ymin');  %------------------(celda de profundidad maxima y minima...lecho)implementado por renzo
% élimination des variables auxiliaires

clear velnan;

