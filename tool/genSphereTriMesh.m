%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           generate TriMesh of Sphere
%           output: TriMesh
%
%           description: TriMesh contact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
format long;
%% Parameters
radius = 10.0;
num = 50;
steps = 1000;
G = 1e-1;
dt = 1.0;
ds = 0.01;
meshFile = 'meshFile.txt';
pointFile = 'pointFile.txt';
%% Points
points = zeros(num,3);
for ii = 1:num
    theta = rand()*pi-pi/2;
    phi = rand()*2*pi;
    x = radius*cos(theta)*cos(phi);
    y = radius*cos(theta)*sin(phi);
    z = radius*sin(theta);
    points(ii,:) = [x y z];
end
%% Uniform
figure (1)
hold on;
for ii = 1:steps
    force = zeros(num,3);
    for jj = 1:num-1
        for kk = jj:num
            dist = points(kk,:) - points(jj,:);
            F = G/((max(norm(dist),ds))^3)*dist;
            force(jj,:) = -F + force(jj,:);
            force(kk,:) =  F + force(kk,:);
        end
    end
    for jj = 1:num
        F = force(jj,:);
        force(jj,:) = F - dot(F,points(jj,:)/norm(points(jj,:)))*points(jj,:)/norm(points(jj,:));
    end
    points = force*dt + points;
    if rem(ii,10)==0 
        plot3(points(:,1),points(:,2),points(:,3),'o') 
    end;
    axis equal;
end
%% Delaunay
DT = delaunayTriangulation(points);
[T,Xb] = freeBoundary(DT);
TR = triangulation(T,Xb);
trisurf(T,Xb(:,1),Xb(:,2),Xb(:,3), ...
     'FaceColor', 'cyan', 'faceAlpha', 0.8);
axis equal;
hold on;
% Calulate the incenters and face normals.
P = incenter(TR);
fn = faceNormal(TR);
% Display the normal vectors on the surface.
quiver3(P(:,1),P(:,2),P(:,3), ...
     fn(:,1),fn(:,2),fn(:,3),0.5, 'color','r');
hold off;
%% Write data
triPoints = TR.Points;
triMesh = TR.ConnectivityList - 1;
fid = fopen(pointFile,'wt');
for ii = 1:size(triPoints,1)
    fprintf(fid,'%12.4e %12.4e %12.4e\n',triPoints(ii,:));
end
fid = fopen(meshFile,'wt');
for ii = 1:size(triMesh,1)
    fprintf(fid,'%4i %4i %4i %4i\n',3,triMesh(ii,:));
end
fclose('all');
