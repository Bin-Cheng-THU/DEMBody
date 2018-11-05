%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           scale TriMesh
%           output: TriMesh
%
%           description: X轴 -> Y轴 -> Z轴， 由负到正
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
format long;
%% Input
RadiusO = 220;
RadiusS = 260;
inputname = 'Bennu.vtk';
%% Data
[points,faces] = loadVTK(inputname);
%% Write
filename = [[inputname(1:end-4),num2str(RadiusS)],'.vtk'];
writeVTK(points/RadiusO*RadiusS,faces,filename);