%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) TU Delft All Rights Reserved
% Written by Abbas Alloul, Giorgio Gardella and Ali Moradvandi
% For any correspondence: a.moradvandi@tudelft.nl

%% Introduction of code (purpose)
% ODE implementation of the mechanstic model of purple bactria (PBM) 
% in an open raceway reactor.
%
% The model is discussed in the paper entitiled "A novel mechanistic 
% modelling approach for microbial selection dynamics:
% towards improved design and control of raceway reactors
% for purple bacteria".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Structure of the code
% Invoke parameters and variables
% Run the model
% Plot the results

clear all
clc

%% Invoke parameters and variables
%% Input data

Days        = 30;                       % Simulation duration (day)                  
Timesteps   = 24;                       % Timestep for each day (e.g. 24 means every hour)
swv         = 1;                        % Switching
V           = 100;                      % Reactor volume (L)
Indata      = Influent';                % Influent characteristic
steps       = 0:1:Days*Timesteps-1;     % 
intensity   = 54;                       % Light intensisty (W/m2)

initial     = [0, 0, 3000, 3.57, 0, 314,  182 ,0, 80, 10, 0, 80,80, 0, 0, V];   % Initial conditions 
            % [O2, SS, SVFA, SIC, SH2, SIN, SIP, SI, XPB_ph, XPB_ch, XPB_an, XAHB, XAN, XS, XI, V]


%% REACTOR GEOMETRY
HRT         = 2.0;                      % day
SRT         = 2.0;                      % day
A           = 0.5;                      % Area, m^2 (h = 0.20m -->20 cm )
h           = V/A/1000;                 % Hight, m
fHS         = HRT/SRT;                  % HRT/SRT ratio defines the fraction of removed particles

%% Schedule flow rate, light and paddlewheel
% Flow rate is defiended for a sequencing batch conditions. It can be
% changed for continouse condition.
Qout_h          = zeros(24,1);          
Qin_h           = zeros(24,1);
Qout_h(24)      = V/HRT;                % Start effluent extraction SBR
Qin_h(24)       = V/HRT;                % Start influent filling SBR

% Paddlewheel is considered working only during days. It can modefied for
% different operation conditions.
sw_h            = ones (24,1); 
sw_h(12:end)    = 0;                    % If "1" paddlewheel on during night

% Time light is considred 12 hours. It can be adjusted for different
% operation conditions.
timelight       = zeros(24,1);
timelight(1:12) = intensity;

% Data organization
Qin             = zeros(length(steps),1);
Qout            = zeros(length(steps),1);
sw              = zeros(length(steps),1);
light_T         = zeros(length(steps),1);

for i = 0:(Days-1)
    Qin((1+(24*i)):(24*(i+1)))      = Qin_h; 
    Qout((1+(24*i)):(24*(i+1)))     = Qout_h; 
    sw((1+(24*i)):(24*(i+1)))       = sw_h;
    light_T((1+(24*i)):(24*(i+1)))  = timelight;
end

%% Run model
% ODE input variables
Input   = [light_T,Qin,Qout,sw];      
options = odeset('NonNegative',1:14);           
[t, y]  = ode15s(@(t,y) PBM(t,y,Indata,Input,fHS,h), steps, initial, options);

%% Plot results 

rate    = zeros(length(steps),10);     % Rates of growth 
I       = zeros(length(steps),11);     % Inhibition parameters
SE      = zeros(length(steps),1);      % Light in water
Yield   = zeros(length(steps),3);      % Flow in/out/volume
Qi      = zeros(length(steps),3);      % Biomass yields 
for i=1:length(y)
    [yy,Qii,Yi]                 = PBM(steps(i),y(i,:),steps,Input,fHS,h);                 
    Qi(i,1:3)                   = Qii;       
    Yield(i,1:3)                = Yi;     
 end


CODeff      = 100 - (y(:,2) + y(:,3))./(Input(:,2) + Input(:,3))*100;                                       % SCOD removal efficiency
CODrr       = ((Input(:,2) + Input(:,3)) - (y(:,2) + y(:,3)))/HRT/1000;                                     % COD removed g/L/d
CODprod     = (y(:,9) + y(:,10) + y(:,11) + y(:,12) + y(:,13) + y(:,14) + y(:,15))/SRT;
Protprod    = ((y(:,9) + y(:,10)+ y(:,11))*.086 + (y(:,12) + y(:,13))*.07)*6.25/SRT;
XPBp        = (y(:,9) + y(:,10) + y(:,11))./(y(:,9) + y(:,10) + y(:,11) + y(:,12) + y(:,13));               % Percentage of XPB
XAHBp       = y(:,12)./(y(:,9) + y(:,10) + y(:,11) + y(:,12) + y(:,13));                                    % Percentage of XAHB
XANp        = y(:,13)./(y(:,9) + y(:,10) + y(:,11) + y(:,12) + y(:,13));                                    % Percentage of XAN
CODr        = [CODeff, CODrr];
Xp          = [XPBp,XAHBp,XANp];

plot(steps, y, steps, Xp, steps, CODr, steps, Yield, 'linewidth',1.3);
grid on
plotbrowser('on');
% Set legend
l = legend('$O_2 (mgO2/L)$','$S_S (mgCOD/L)$','$S_{VFA} (mgCOD/L)$','$S_IC (mgCOD/L)$',...
    '$S_{H2}(mgCOD/L)$','$S_{IN} (mgN/L)$','$S_{IP} (mgP/L)$','$S_I (mgCOD/L)$',...
    '$X_{PB,ph} (mgCOD/L)$','$X_{PB,ch} (mgCOD/L)$','$X_{PB,an} (mgCOD/L)$',...
    '$X_{AHB} (mgCOD/L)$','$X_{AN} (mgCOD/L)$','$X_S (mgCOD/L)$','$X_I (mgCOD/L)$',...
    '$V (L)$',...
    '$X_{PBp} (\%)$','$X_{AHBp} (\%)$','$X_{ANp} (\%)$',...
    '$COD_{eff} (\%)$', '$COD_{rr}$', '$yield_{PPB}$', '$yield_{AHB}$', '$yield_{tot}$');
set(l,'Interpreter','Latex','fontsize',10);
% Set title
t = title('Output purple phototrophic model');
set(t,'Interpreter','Latex','fontsize',12);
% Set axis labels
xlabel('Time [$h$]','Interpreter','Latex','fontsize',12);

