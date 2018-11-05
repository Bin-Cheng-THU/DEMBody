function  writeVTK(points,faces,filename)
%
%specified type: triangle facet polyhydron vtk file
%vectors all represented in body coordinates, default unit: m
disp('Start writting...');
fid = fopen(filename,'w+');
fprintf(fid,'# vtk DataFile Version 3.0\r\n');
fprintf(fid,'%s\r\n',filename);
fprintf(fid,'ASCII\r\n');
fprintf(fid,'DATASET POLYDATA\r\n');
fprintf(fid,'POLYGONS %d %d\r\n',length(faces),length(faces)*4);
for ii = 1:length(faces)
    fprintf(fid,'%d %d %d %d\r\n',3,faces(ii,:));
end
fprintf(fid,'POINTS %d float\r\n',length(points));
for ii = 1:length(points)
    fprintf(fid,'%f %f %f\r\n',points(ii,:));
end
fclose(fid);
disp('Writting completed.')
end

