%% Engineering optimization 
% Governing parameters
% Martin Janssens
% Rens Liebrand

%% General 
L=1.65;                         %m
P_skier_max=780;                %N
P_skier=400;                    %N
G=9.81;                         %m/s
nu_snow=0.05;                   %Pa/s

%% Aerodynamic
rho_air=1.225;                  %kg/m^3
Uwind=4;                        %m/s
Uskier=10;                      %m/s
Uinduced=2;                     %m/s
Uinf=Uwind+Uskier+Uinduced;     %m/s
nu=1.81206e-5;                  %Pa/s
mu=1.5e-5;                      %m^2/s

%% Structural
rho_carbon=1638;                %kg/m^3
E=200e9;                        %Pa
K=1;                            % - 