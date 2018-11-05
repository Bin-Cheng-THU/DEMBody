%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           generate Gravity TriMesh
%           output: TriMesh
%
%           description: X轴 -> Y轴 -> Z轴， 由负到正
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
format long;
%% Input
Dx = [-400 400];
Dy = [-400 400];
Dz = [-400 400];
grid = 10;
filename = 'particles.txt';
isParaview = 1;
%% Data
Nx = (Dx(2)-Dx(1))/grid+1;
Ny = (Dy(2)-Dy(1))/grid+1;
Nz = (Dz(2)-Dz(1))/grid+1;

points = zeros(Nx*Ny*Nz,3);
index = 1;
for ii = 1:Nz
    point(3) = Dz(1) + (ii-1)*grid;
    for jj = 1:Ny
        point(2) = Dy(1) + (jj-1)*grid;
        for kk = 1:Nx
            point(1) = Dx(1) + (kk-1)*grid;
            points(index,:) = point;
            index = index +1;
        end
    end
end
%% Write Data
fid = fopen(filename,'wt');
fprintf(fid,'%12i\n',length(points));
for ii = 1:size(points,1)
    fprintf(fid,'%12.4e %12.4e %12.4e\n',points(ii,:));
end
fclose(fid);
%% Paraview File
if (isParaview)
    paraname = [filename(1:end-4),'Para.csv'];
    fid = fopen(paraname,'wt');
    fprintf(fid,'X, Y, Z\n');
    for ii = 1:size(points,1)
        fprintf(fid,'%12.4e, %12.4e, %12.4e\n',points(ii,:));
    end
    fclose(fid);
end










