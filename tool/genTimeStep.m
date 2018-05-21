%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           generate Time step suitable for corresponding material
%           input: material properties
%           output: time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear
format long;
%% Input parameters
E = 0.25e8; % Young's Modules
Rho = 2.15; % Density
nu = 0.20; % Poisson's ratio
v = 100.0; % maximum velocity
r = 0.5; % average radius

%% Estimate time step
T = 5.84*(Rho*(1-nu^2)/E)^(2/5)*r*v^(-1/5);
dT = T/20;
disp(dT);

%% Compare to Linear Model
Kn = 2.0*E*sqrt(r/2)/(3*(1-nu^2));
disp(Kn)

MaxAcc = 2.0*E*sqrt(r)/(3*(1-nu^2))*(r)^(3/2)/(4/3*pi*r^3*Rho);
disp(MaxAcc)