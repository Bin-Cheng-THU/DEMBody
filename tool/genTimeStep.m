%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           generate Time step suitable for corresponding material
%           input: material properties
%           output: time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear
format long;
%% Input parameters
global E Rho nu v r  m
E = 1e4; % Young's Modules
Rho = 0.4*1000; % Density
nu = 0.30; % Poisson's ratio
v = 0.01; % maximum velocity
r = 5.0; % average radius
m = 4/3*pi*r^3*Rho; % particles mass
%% Estimate time step
T = 5.84*(Rho*(1-nu^2)/E)^(2/5)*r*v^(-1/5);
dT = T/20;
disp(dT);
%% Compare to Linear Model
Kn = 2.0*E*sqrt(r/2)/(3*(1-nu^2));
disp(Kn)

MaxAcc = 2.0*E*sqrt(r*r)/(3*(1-nu^2))*(r)^(3/2)/(4/3*pi*r^3*Rho);
disp(MaxAcc)
%% Numerical simulation
[Time,Pos]=ode45('contact_model',[0 T*2],[0 v]);
plotyy(Time,Pos(:,1),Time,Pos(:,2));