%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           convert VTK mesh files to DEMBody mesh files
%           input: mesh file; point file
%           output: DEMBody file
%
%           description: TriMesh contact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
format long;
%% Parameters
%---file name
filename = 'Bennu.vtk';
outputfile = strcat(filename(1:end-3),'mesh');
%% Load Data
[points,faces] = loadVTK(filename);
TR = triangulation(faces+1,points);
%% Triangle analyse
figure (1)
trisurf(TR.ConnectivityList,points(:,1),points(:,2),points(:,3), ...
     'FaceColor', 'cyan', 'faceAlpha', 0.8);
axis equal;
hold on;
% Calulate the incenters and face normals.
P = incenter(TR);
fn = faceNormal(TR);
% Display the normal vectors on the surface.
quiver3(P(:,1),P(:,2),P(:,3), ...
     fn(:,1),fn(:,2),fn(:,3),0.5, 'color','r');
hold on;
%% Generate data
num = size(faces,1);
trimeshPoint = zeros(num,3);
trimeshVectorN = zeros(num,3);
trimeshVectorTx = zeros(num,3);
trimeshVectorTy = zeros(num,3);
for ii = 1:num
    pointA = points(faces(ii,1)+1,:);
    pointB = points(faces(ii,2)+1,:);
    pointC = points(faces(ii,3)+1,:);
    AB = pointB - pointA;
    AC = pointC - pointA;
    % Check the order
    if dot(cross(AB,AC),fn(ii,:))<0.0
        disp('Error')
        pointB = points(faces(ii,3)+1,:);
        pointC = points(faces(ii,2)+1,:);
        AB = pointB - pointA;
        AC = pointC - pointA;
    end
    trimeshPoint(ii,:) = pointA;
    trimeshVectorTx(ii,:) = AB;
    trimeshVectorTy(ii,:) = AC;
    vectorN = cross(AB,AC);
    vectorN = vectorN/norm(vectorN);
    trimeshVectorN(ii,:) = vectorN;
%     L1 = dot(AC,AC)*dot(AB,AB) - dot(AB,AC)*dot(AC,AB);
%     L2 = dot(AC,AC);
%     L3 = dot(AC,AB);
%     L4 = dot(AB,AB);
%     trimeshLength(ii,:) = [L1,L2,L3,L4];
end

%% Write data
fid = fopen(outputfile,'wt');
for ii = 1:num
    fprintf(fid,'%12.4e %12.4e %12.4e',trimeshPoint(ii,:));
    fprintf(fid,'%12.4e %12.4e %12.4e',trimeshVectorTx(ii,:));
    fprintf(fid,'%12.4e %12.4e %12.4e\n',trimeshVectorTy(ii,:));
end
fclose('all');
    
