%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           generate Mesh Parameters for DEMBody
%           input: particle size; boundary size
%           output: DEMBody Parallel Lattice Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
format long;
%% Input Parameters
dx = 0.75; %2.5Rmax
x = 24.0; %Boundary size
y = 24.0; %Boundary size
z = 24.0;
NMAX = 500000;

%% Mesh Grid
Nx = x/dx*2;
Ny = y/dx*2;
num = Nx*Ny;

casenum(1) = 1;
casenum(2) = -Nx-1;
casenum(3) = -Nx;
casenum(4) = -Nx+1;
casenum(5) = -num;
casenum(6) = -num+1;
casenum(7) = -num+Nx-1;
casenum(8) = -num+Nx;
casenum(9) = -num+Nx+1;
casenum(10) = -num-1;
casenum(11) = -num-Nx+1;
casenum(12) = -num-Nx;
casenum(13) = -num-Nx-1;
fprintf(1,'%8d,%8d,%8d,%8d,%8d,%8d,%8d,%8d,%8d,%8d,%8d,%8d,%8d\n',abs(casenum(1:13)));
disp(Nx);
disp(num);

%% Parallel Lattice