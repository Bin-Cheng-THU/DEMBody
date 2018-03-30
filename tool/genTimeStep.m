%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           generate Time step suitable for corresponding material
%           input: material properties
%           output: time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear
format long;
%% Input parameters
E = 0.25e9; % Young's Modules
Rho = 2.15; % Density
nu = 0.20; % Poisson's ratio
v = 10; % maximum velocity
r = 0.5; % average radius
%% Estimate time step
T = 5.84*(Rho*(1-nu^2)/E)^(2/5)*r*v^(-1/5);
dT = T/50;
disp(dT);
