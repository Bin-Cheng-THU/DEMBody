function  writeBT(points,faces,filename)
%
%specified type: triangle facet polyhydron vtk file
%vectors all represented in body coordinates, default unit: m
disp('Start writting...');
fid = fopen(filename,'w+');
fprintf(fid,'%d %d\r\n',length(points),length(faces));
for ii = 1:length(points)
    fprintf(fid,'%f %f %f\r\n',points(ii,:));
end
for ii = 1:length(faces)
    fprintf(fid,'%d %d %d\r\n',faces(ii,:));
end
fclose(fid);
disp('Writting completed.')
end

