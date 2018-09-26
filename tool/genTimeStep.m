%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           generate Time step suitable for corresponding material
%           input: material properties
%           output: time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear
format long;
%% Input parameters
global E Rho nu v r  m
E = 0.1e7; % Young's Modules
Rho = 2.5; % Density
nu = 0.29; % Poisson's ratio
v = 10.0; % maximum velocity
r = 0.2; % average radius
m = 4/3*pi*r^3*Rho; % particles mass
%% Estimate time step
T = 5.84*(Rho*(1-nu^2)/E)^(2/5)*r*v^(-1/5);
dT = T/20;
disp(dT);
%% Compare to Linear Model
Kn = 2.0*E*sqrt(r/2)/(3*(1-nu^2));
disp(Kn)

MaxAcc = 2.0*E*sqrt(r)/(3*(1-nu^2))*(0.2*r)^(3/2)/(4/3*pi*r^3*Rho);
disp(MaxAcc)
%% Numerical simulation
[Time,Pos]=ode45('contact_model',[0 0.005],[0 v]);
plotyy(Time,Pos(:,1),Time,Pos(:,2));