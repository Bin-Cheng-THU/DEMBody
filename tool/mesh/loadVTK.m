function [points,faces] = loadVTK(filename)
%
%specified type: triangle facet polyhydron vtk file
%vectors all represented in body coordinates, default unit: m
disp('Start preloading...');
fid = fopen(filename);
tline=fgetl(fid);
tline=fgetl(fid);
tline=fgetl(fid);
tline=fgetl(fid);
tmp = fscanf(fid,'%*s %d %d',[1,2]);
nFaces = tmp(1);
faces = fscanf(fid,'%d %d %d %d',[4,nFaces]);
tmp = fscanf(fid,'%*s %d %*s',[1,1]);
nPoints = tmp;
tline=fgetl(fid);
points = fscanf(fid,'%f %f %f',[3,nPoints]);
fclose(fid);
disp('Preloading completed.')
end

