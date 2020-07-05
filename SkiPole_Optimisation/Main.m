%% Engineering optimization 
% Main
% Martin Janssens
% Rens Liebrand
clear all
close all
warning off

% Import all system constants
Params;

% Import database
F=import_database();

% Parameters in middle of range
% In CENTIMETERS
% 0.5 < c < 4
% 0.2c < h < 0.45c           ->   0.05 < h < 1.35
% 0.3c < r < 0.5c           ->   0.15 < r < 2
% 0.08 < t < 0.1h=0.2*0.45*c  ->   0.08 < t < 0.48

% Set avg parameters
t=3e-3;
r=0.4; % r   *c
h=0.25; % h   *c
c=0.012;

%% Running Different part of the exercise
% run Prelim_analysis
% run Problem_Investigation
% run Simplified_Problem
% run Optimisation_Fitted_Response
% run Own_Optimisation_SD
% run Optimisation_Fitted_Response_Gradient_Projection
run Full_Optimisation_Smooth















