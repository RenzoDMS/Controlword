function  ROUSE2
close all
clc
% Insert below the file you wish to read
data = load('C:\Users\lapa-8\Documents\Philemon\MATLAB\Sedim_CONC\ConcentrationTEMPLATE2.txt'); 
depth = data(:,1);% depth
C_obs = data(:,4);% observed concentrations
maxdepth = 25.50;% max depth CHANGE THIS VALUE
h=maxdepth;
a=0.05*maxdepth;
z0=maxdepth - depth;

%Initiate search parameters
%Put an approximate value of W and Cz0
W = 1e-4;
Cz0 = 10; %[mg/l]
X0 = [Cz0,W];

%Function to optimise the parameters using an objective function defined
%below
[X OF] = fminsearch(@Objfunction, X0, [], C_obs, h, a, z0);
%Extract optimised parameters
Cz0opt = X(1);
Wopt = X(2);

%Plot observed and modelled concentrations
z0_plot =0:0.1:h;
Ccalc_plot = rouse(Cz0opt,Wopt,h,a,z0_plot);
plot(C_obs, z0, 'x', 'Color', 'k');
hold on
plot(Ccalc_plot, z0_plot, 'r');
%set(gca,'YTick',linspace(-maxdepth,0,5));
hold on
plot(0:max(C_obs)+10, repmat(h, numel(0:max(C_obs)+10),1), 'b-', 'LineWidth', 1.2);
hold on
plot(max(C_obs),h+0.5, 'v', 'MarkerFaceColor', 'b', 'MarkerSize',8)
xlabel('Concentrations [mg/L]')
ylabel('Elevation above the river bed [m]  ')
title('Observed and modelled concentrations of fines using the Rouse profile')
hold off
%Display optimized parameters
fprintf('\nOF = %6.3f',OF);
fprintf('\nOptimised W = %f ',Wopt);
fprintf('\nOptimised Cz0    = %f [mg/L] \n',Cz0opt);

%Objective function definition
function OF= Objfunction(X, C_obs, h, a, z0)
Cz0 = X(1);
W = X(2);
C_calc = rouse(Cz0,W,h,a,z0);
OF1=mean(abs(C_calc-C_obs));
OF2=sqrt(mean((C_calc-C_obs).^2));
OF3=mean(abs(log(C_calc)-log(C_obs)));
OF4=mean(abs((C_calc-C_obs)./C_obs));
%SELECT HERE THE OBJECTIVE FUNCTION YOU WANT TO USE.
%Choose between:
%OF1 = mean relative error
%OF2 = RMSE
%OF3 = log of error
%OF4 = mean relative error
OF=OF2; % Change here according to objective function you want to use
end

%Rouse Function DO NOT EDIT

    function Cf = rouse(Cz0,W,h,a,z0)
        Cf = Cz0.*(((h-z0)./z0).*(a./(h-a))).^W;
    end
end


