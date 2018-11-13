%%----------------------------------------------------------------------
% r�cup�ration de la matrice des vitesses � partir du fichier initial
% consitut� de la matrice des profondeurs + matrice des vitesses agr�g�es
% l'une apr�s l'autre.
%----------------------------------------------------------------------
clear all;
close all;
clc;
parametros=[0 0 0 0 0];

% demande � l'�cran de verticales � analyser
  aforos = input('Ingrese numero de aforos: ');
  
% Empezamos a analizar las tres verticales de medidas realizadas
% k=es la cantidad de verticales analizadas

k=1;
for k=1:aforos
% demande � l'�cran de verticales � analyser
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
%calcul du nombre de lignes � rajouter en surface
hauteur_cellule = prof_temp(2,1)-prof_temp(1,1);
nb_cel_ajout = floor(prof_temp(1)/hauteur_cellule);
noe=A.Sup.noe; % number of ensembles
[nbligtemp nbcol] = size(prof_temp);
nblig = nb_cel_ajout + nbligtemp;
%cr�ation des matrices avec les lignes � rajouter en plus
prof = ones(nblig, nbcol);
velocity = ones(nblig, nbcol);

%calcul des profondeurs sur les cellules sup�rieures
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
% �limination des variables auxiliaires
%clear data;
%clear aux_vel;
%clear pos_max_deep;
%clear prof_temp;
%clear velocity_temp;

%%
%--------------------------------------------------------------------
% cr�ation du vecteur de profondeur max o� une vitesse a �t� enregistr�e
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
nbcol = t(:,2); % nombre de verticales � traiter 
nblig = t(:,1); % nombre de cellules max ADCP par ensemble
%idp=[]
pmax=ones(1,nbcol);
d=1;              
% cr�ation du vecteur de profondeur max avec valeur de vitesse
velnan=isnan(velocity);
  for i=d:nbcol,
       pmax(i)=find(velnan(:,i)==0,1,'last');
    % eval(['pmax(' num2str(i) ')=id;']);
  end

% cr�ation du vecteur de profondeur min sans valeur de vitesse
pmin=ones(1,nbcol);
  for i=d:nbcol,            
     pmin(i)=find(velnan(:,i)==1,1,'first');
   % eval(['pmin(' num2str(i) ')=id;']);
  end

Pp=[pmax;pmin];

% �limination des variables auxiliaires
clear velnan;

%------------------------------------------------------------------------
% detection des trous de donn�es
% c.a.d. valeurs de zeros entre des valeurs de vitesses non nulles
%------------------------------------------------------------------------
if (max((pmax+1)- pmin)>0)  
    disp('il y a des trous dans les donn�es de vitesse');
    verticale=find(((pmax+1)-pmin));
    disp('aux verticales');
    disp(verticale);
    
    % nettoyage des trous: moyenne des deux voisins sur les colonnes
    % adjacentes
    dmax=length(verticale);
    for i=d:dmax,
        
        if verticale(i)==1 
            res=velocity (pmin(verticale(i)), verticale(i)+1);
            
        else if verticale(i)==nbcol 
                res=velocity (pmin(verticale(i)), verticale(i)-1);
            else res = (velocity (pmin(verticale(i)), verticale(i)+1)+ velocity (pmin(verticale(i)), verticale(i)-1))/2;    
            end
        end
        velocity (pmin(verticale(i)), verticale(i)) = res;
    end
else disp('pas de trou de donn�es observ�es');
end


% �limination des variables auxiliaires
clear verticale;
clear res;

% --------------------------------------------------------------------------
% construction d'une matrice de vitesses avec une r�f�rence � partie du
% fond
%---------------------------------------------------------------------------
% construction d'une matrice de profondeur 
zblanche = 0.06; % epaisseur de la zone blanche ADCP

pfond = prof(pmax,1)+0.3; %pfond = prof(pmax,1)/(1-zblanche); se modifrico por renzo por altura del bin.

% creation de la matrice flip
velocity_flip = NaN(size(velocity));
% trouver pour chaque verticale les vitesses non nulles entre la surface et
% le fond et 

d=1;
for i=d:nbcol,
    id=find(prof(:,i)<= pfond(i)); %find(round(prof(:,i))<= round(pfond(i))) se modific�, por renzo
    vtemp = velocity(id,i);
    vtemp2 = flipud(vtemp);
    velocity_flip(1:length(vtemp2),i) = vtemp2;
end

num_figure = num_figure+1;
figure(num_figure);
%imagesc(velocity_flip);
imagesc(distance_temp, prof(:,1),velocity_flip);
title(namedate,'FontName','Arial','FontSize', 14);
ylabel('depth (m)','FontName', 'Arial','FontSize', 14);
xlabel('distance (m)','FontName', 'Arial','FontSize', 14);
colorbar

velocity_flip(1:3,:)=NaN;%utilizar solo quando tenga errores en las 
%velocidades de inicio
% �limination des variables auxiliaires
clear zblanche;
clear vtemp;
clear vtemp2;
%%
%----------------------------------------------------------------------
% calcul de vitesse moyenne par verticale choisie par l'utilisateur
% en augmentant la fen�tre de moyennisation de 1 en 1
%----------------------------------------------------------------------



% demande � l'�cran de verticales � analyser
  dist = input('distance de la verticale � moyenner [dist1 dist2 ... distN]: ');

 % calcul de verticales de vitesses moyennes
  vitmoy = ones(nblig ,length(dist),nbcol);
  lmax = ones(length(dist));
  incrementmax = 200; % taille maximale de la fen�tre de  moyennisation

  for i=1:length(dist)
    % recherche de la position de la verticale
      mdif= min(abs(distance(:)- dist(i)));
      posdist = min(find(mdif==abs(distance(:)- dist(i)))); % -------------- posdist=min(find(abs(round(distance(:))- dist(i))<= 1.5)); se ajust�, por Renzo.
        
    % longueur maximale de moyennisation
      lmax(i)= min([(posdist-1) (nbcol-posdist) 200]);
    
      for inc=1:lmax
        % isoler la zone de verticales autour de la verticale choisie
          zonev = velocity_flip(:,posdist-inc:posdist+inc);
          vitmoy(:,i,inc) = nanmean(zonev,2);
      end
    end


%
%------------------------------------------------------------------------
% impression des verticales de vitesse moyennes pour une taille de fen�tre
% fixe 
%------------------------------------------------------------------------
% definition des verticales moyennes � afficher
  num_figure = num_figure+1;
  figure(num_figure);

  incplot = [1 round(lmax(i)/10) round(lmax(i)/2) lmax(i)];   %se puede modificar la cantidad de verticales para el grafico aqui
  nVelprome=incplot
  cor={'ko-';'co-';'bo-';'go-';'m-o';'y-o'};

  coefhauteur= 0.8 ;

  for i=1:length(dist)
      if length(dist) ~= 1 
          subplot(2,length(dist)/2,i); 
      end;
    
      hold on;
    
      for inc = 1:length(incplot)
          vitmoytemp = vitmoy(:,i,incplot(inc));       
          id = find(isnan(vitmoytemp)==0,1 ,'last');
          prof_coehauteur=round(prof(id)*coefhauteur);
          id_prof=find(prof(:,1)<=prof_coehauteur);
        
          plot(vitmoy(1:id_prof(end),i,incplot(inc)),prof(1:id_prof(end)),cor{inc});%grafico
        
        %com profundidade total
       % plot(vitmoy(1:id,i,incplot(inc)),prof(1:id_prof),cor{inc});%grafico com profundida en funsi�n de CoefHauteur
        %xlim([0.5 1.2*max(vitmoy(:,i, incplot(inc)))]);
      end 
    
      title('verticale','FontName', 'Arial','FontSize', 12);
      xlabel('velocity (m/s)','FontName', 'Arial','FontSize', 12);
      ylabel('depth (m)','FontName', 'Arial','FontSize', 12);
      legend(sprintf('%d',incplot(1)),sprintf('%d',incplot(2)),sprintf('%d',incplot(3)), sprintf('%d',incplot(4)));
  end
  hold off
% �limination des variables auxiliaires
%clear vitmoytemp;

%------------------------------------------------------------------------
% calcul des parametres de regression d'une courbe log sur une verticale 
% de vitesse moyennes sur une hauteur d'eau de 0.8*H et affichage de
% l'�volution de ces param�tres en fonction de la taille de la fen�tre de
% moyennisation.
%------------------------------------------------------------------------
%a = zeros(max(lmax),length(dist)); % tableau de valeur du coefficient du log
%b = zeros(max(lmax),length(dist)); % tableau de valeur de la constante du log
  coefhauteur = 0.8 ; % pourcentage de hauteur sur laquelle la r�gression est effectu�e
  for i=1:length(dist)
      for inc=1:lmax(i)
          vitmoy_temp = vitmoy(:,i,inc);       
          idmax = find(isnan(vitmoy_temp)==0,1 ,'last');
          hauteurmax = prof(idmax);
          hauteurfit = coefhauteur*hauteurmax;
          id_hauteurfit = min(find(abs(prof(:,i)- hauteurfit)<= hauteur_cellule*1.5));
          x = prof(1:id_hauteurfit,1);
          v = vitmoy(1:id_hauteurfit,i,inc);
        
        % function [slope, intercept,MSE, R2, S] = logfit(x,y,varargin) ;
        % v=aln(x)+b=(a/log10(e))*log10(x)+b;
      
          [slope, intercept, MSE, R2, S, extra]= logfit(x,v,'logx');
        
          a(inc,i) = slope; %pendiente para logaritmo base 10
          b(inc,i) = intercept;% intercepto
          r2(inc,i)= R2; %grado de correlacion
      end
  end
  
%-------------------------------------------------------------------------------------
%hallando la profundidad maxima en posdist y espesor de fondo;
  hmax=pfond(posdist)
  Z0=0.05*hmax
%identificando la mejor regresi�n.
  disp('mejor regresi�n v= ap*log(x)+ b= (ap*log10(e))*ln(x)+ b ');
  nperfsum=find(r2==max(r2))
  ap=a(nperfsum)
  bp=b(nperfsum)
  r2=max(r2)
  vitmoy_temp = vitmoy(:,i,nperfsum);       
          idmax = find(isnan(vitmoy_temp)==0,1 ,'last');
          hauteurmax = prof(idmax);
          hauteurfit = coefhauteur*hauteurmax;
          id_hauteurfit = min(find(abs(prof(:,i)- hauteurfit)<= hauteur_cellule*1.5));
          xp = prof(1:id_hauteurfit,1);
          vm= vitmoy(1:id_hauteurfit,i,nperfsum); %velocidad muestreada
%vamos obtener vector velocidad con la curva ajustada
  for j=1:length(xp)
      veloc=ap*log10(xp(j))+ bp; %funcion velocidad ajustada
      vp(j,1)=veloc;
  end
  %muestra=[muestra xp vm];
  parametro=[ap bp r2 hmax Z0];
  parametros=[parametros; parametro];
  
  kk=dist;
  %xlswrite('velocity_temp.xlsx', muestra ,sprintf('xp vm V%d ',kk ));  %------------------(vector a b r2 )implementado por renzo
  xlswrite('velocity_temp.xlsx', parametros , sprintf('ap bp r2 V%d ',kk ));  %------------------(vector a b r2 )implementado por renzo
%grafica de la velocidad vs altura para el mejor ajuste
  num_figure = num_figure+1;
  figure(num_figure);
  if xp(1)<=0.0001;
      xp(1)=NaN;
      vp(1)=NaN;
  end
  plot(vitmoy(1:id_hauteurfit,1,nperfsum), prof(1:id_hauteurfit),'ko-',vp(1:size(xp)) , xp(1:size(xp)) ,'co-');%grafico
  aa=ap*log10(2.718281828459045235360);
  aa=round(aa,4);
  bb=round(bp,4);
      title('Compare','FontName', 'Arial','FontSize', 12);
      xlabel('velocity (m/s)','FontName', 'Arial','FontSize', 12);
      ylabel('depth (m)','FontName', 'Arial','FontSize', 12);
      legend('o sample', sprintf('o fit_v: \n v(x)=%f ln(x)+ %f',aa,bb));
%--------------------------------------------------------------------------------------

%for inc = 1:length(incplot)
%pp=60;
%x2 = [1:1:55];
%f = a(incplot(inc),i)*log(x2)+ b(incplot(inc),i);
%plot(f,x2,'r-');
%end
%hold off;

% affichage de a vs taille de la fen�tre de moyennisation
  num_figure=num_figure+1;
  figure(num_figure);
  id = find(a(:,1)); % trouve la derni�re valeur de a non nulle
  ax = [1:1:length(id)]; % construction du tableau des num�ro d'incr�mentation
%plot(ax,a(1:max(id)),'r-');
  plot(ax,a,'-');
  title('coefficient(a) vs acumulative mean cell','FontName','Arial','FontSize', 14);
  ylabel('a (m/s)','FontName', 'Arial','FontSize', 14);
  xlabel('No. celulas incrementadas','FontName', 'Arial','FontSize', 14);
  ylabel('a (m/s)','FontName', 'Arial','FontSize', 14);
  legend('a_v1','a_v2','a_v3','a_v4');

% affichage de b vs taille de la fen�tre de moyennisation
  num_figure=num_figure+1;
  figure(num_figure);
  id = find(b(:,1)); % trouve la derni�re valeur de a non nulle
  bx = [1:1:length(id)]; % construction du tableau des num�ro d'incr�mentation
%plot(bx,b(1:max(id)),'b-');
  plot(bx,b,'-');
  title('coefficient(b) vs acumulative mean cell','FontName','Arial','FontSize', 14);
  ylabel('b (m/s)','FontName', 'Arial','FontSize', 14);
  xlabel('No.celulas incrementadas','FontName', 'Arial','FontSize', 14);
  legend('b_v1','b_v2','b_v3','b_v4');

k=k


end


