%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           convert VTK files to POV-Ray files
%           input: mesh file; point file; bondedWall file
%           output: POV-Ray files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear
format long;
%% Parameters
%---folder path of particle datas
folderParticle = 'C:\Users\chengbin\Desktop\Data\Data\';
%---number of files
Num = 200;
Step = 1;
%---camera settings
bgclr = [0.3, 0.3, 0.3];
angle = 4.5;
location = [-40.0, 0.0, 5.0];
skyvec = [0,0,1];
focus = [0,0.0,5.0];
lightsrc = [-8.0,0.0,50.0];
%% Read
tmp = importdata('BondedWalls.txt');
bondedWallX = tmp(1:Step:end,[2 3 4]);
bondedWallQ = tmp(1:Step:end,[11 12 13 14]);
clear tmp;

tmp = importdata('BondedWallMesh.txt');
bondedWallMesh = tmp(:,[2 3 4]);
clear tmp;

tmp = importdata('BondedWallPoint.txt');
bondedWallPoint = tmp;
clear tmp;
%% Pov-Ray file generating
for I = 1000:(1000+Num-1)
    %% read particle data
    filename = strcat(folderParticle,num2str(I));
    filename = strcat(filename,'X.csv');
    
    tmp = importdata(filename);
    tmp = tmp.data;
    point = tmp(:,1:3);
    r = tmp(:,10);
    q = tmp(:,19:22);
    
    %% attitude matrix of particles
    for J = 1:length(r)
        matrix(J,1,1) = q(J,1)^2-q(J,2)^2-q(J,3)^2+q(J,4)^2;
        matrix(J,1,2) = 2*(q(J,1)*q(J,2)+q(J,3)*q(J,4));
        matrix(J,1,3) = 2*(q(J,1)*q(J,3)-q(J,2)*q(J,4));
        matrix(J,2,1) = 2*(q(J,1)*q(J,2)-q(J,3)*q(J,4));
        matrix(J,2,2) = -q(J,1)^2+q(J,2)^2-q(J,3)^2+q(J,4)^2;
        matrix(J,2,3) = 2*(q(J,2)*q(J,3)+q(J,1)*q(J,4));
        matrix(J,3,1) = 2*(q(J,1)*q(J,3)+q(J,2)*q(J,4));
        matrix(J,3,2) = 2*(q(J,2)*q(J,3)-q(J,1)*q(J,4));
        matrix(J,3,3) = -q(J,1)^2-q(J,2)^2+q(J,3)^2+q(J,4)^2;
    end
    
    %% attitude matrix of bonded walls
    bondedWallMatrix = quat2dcm([bondedWallQ(I-999,4) bondedWallQ(I-999,1) bondedWallQ(I-999,2) bondedWallQ(I-999,3)]);
   
    %% write files
    filename = strcat('Data\Scene',num2str(I));
    filename = strcat(filename,'.pov');

    fid = fopen(filename,'wt');
    fprintf(fid,'#include "colors.inc"\n');
    fprintf(fid,'#include "shapes.inc"\n');
    fprintf(fid,'#include "textures.inc"\n');
    fprintf(fid,'#include "stones.inc"\n');
    fprintf(fid,'#include "woods.inc"\n');
    fprintf(fid,'#include "glass.inc"\n');
    fprintf(fid,'#include "metals.inc"\n');
    %
    fprintf(fid,'light_source {<%12.4e,%12.4e,%12.4e> color 1.5}\n',lightsrc);
    fprintf(fid,'light_source {<%12.4e,%12.4e,%12.4e> color 0.5}\n',location);
    fprintf(fid,'background { color rgb <%12.4e, %12.4e, %12.4e> }\n', bgclr);
    fprintf(fid,'camera {location <%12.4e,%12.4e,%12.4e> \n', location);
    fprintf(fid,'        sky   <%12.4e,%12.4e,%12.4e> \n', skyvec);
    fprintf(fid,'        look_at <%12.4e,%12.4e,%12.4e> \n', focus);
    fprintf(fid,'        }\n'); 
    %
    fprintf(fid,'#declare Texture_W = \n');
    fprintf(fid,'    texture{ pigment{ color White*0.8} \n');
    fprintf(fid,'        normal { bumps 1 scale 0.025} \n');
    fprintf(fid,'        finish { diffuse 0.9 specular 1} \n');
    fprintf(fid,'        }\n');
    fprintf(fid,'#declare Texture_S = \n');
    fprintf(fid,'            texture{ T_Wood2 scale 1 \n');
    fprintf(fid,'normal { agate 0.5 scale 0.25} \n');
    fprintf(fid,'finish { diffuse 0.9 phong 1 } \n');
    fprintf(fid,'        }\n');
    fprintf(fid,'#declare FelbriggSand = texture {\n');
    fprintf(fid,'            pigment {color rgb < 1, 0.9, 0.65>}\n');
    fprintf(fid,'normal {granite 0.2 scale 0.02}\n');
    fprintf(fid,'finish {brilliance 0.1 specular 0.3 ambient 0.05}\n');
    fprintf(fid,'        }\n');    
    fprintf(fid,'#declare TdG_DarkSand = texture {\n');
    fprintf(fid,'            pigment {agate agate_turb 0.3\n');
    fprintf(fid,'color_map { [ 0.0     rgbft <1.0, 0.816993, 0.479133, 0.0, 0.0> ]\n');
    fprintf(fid,'[ 0.049822  rgbft <0.939815, 0.699074, 0.518519, 0.0,0.0>*0.1 ]\n');
    fprintf(fid,'[ 0.096085  rgbft <0.991333, 0.730905, 0.383667, 0.0,0.0>*0.1 ]\n');
    fprintf(fid,'[ 0.149466  rgbft <0.990941, 0.863911, 0.446526, 0.0,0.0>*0.1 ]\n');
    fprintf(fid,'[ 0.41637  rgbft <0.925926, 0.739198, 0.407407, 0.0,0.0>*0.1 ]\n');
    fprintf(fid,'[ 1.0     rgbft <1.0, 0.971785, 0.729133, 0.0, 0.0> ]}\n');
    fprintf(fid,'turbulence 5.0 lambda 6.0 frequency 6.0 ramp_wave warp {turbulence <5.0, 2.0, 1.0>}}\n');
    fprintf(fid,'normal{granite , 0.2  warp {turbulence <1.0, 5.0, 2.0>}}\n');
    fprintf(fid,'finish {ambient 0.01 specular 0.1}\n');
    fprintf(fid,'        }\n');      
    fprintf(fid,'#declare TdG_LightSand  = texture {\n');
    fprintf(fid,'            pigment {agate agate_turb 0.3\n');
    fprintf(fid,'color_map {[ 0.0     rgbft <1.0, 0.816993, 0.479133, 0.0, 0.0> ]\n');
    fprintf(fid,'[ 0.049822  rgbft <0.939815, 0.699074, 0.518519, 0.0, 0.0> ]\n');
    fprintf(fid,'[ 0.096085  rgbft <0.991333, 0.730905, 0.383667, 0.0, 0.0> ]\n');
    fprintf(fid,'[ 0.149466  rgbft <0.990941, 0.863911, 0.446526, 0.0, 0.0> ]\n');
    fprintf(fid,'[ 0.41637  rgbft <0.925926, 0.739198, 0.407407, 0.0, 0.0> ]\n');
    fprintf(fid,'[ 1.0     rgbft <1.0, 0.971785, 0.729133, 0.0, 0.0> ]}\n');
    fprintf(fid,'turbulence 5.0 lambda 6.0 frequency 6.0 ramp_wave warp {turbulence <5.0, 2.0, 1.0>}}\n');
    fprintf(fid,'normal{granite , 0.2  warp {turbulence <1.0, 5.0, 2.0>}}\n');
    fprintf(fid,'finish {ambient 0.01 specular 0.1}\n');
    fprintf(fid,'        }\n');
    %
    for J = 1:length(r)
        %-------------
        fprintf(fid,'sphere ');
        fprintf(fid,'{<%25.16e,%25.16e,%25.16e>, %25.16e \n',[0 0 0],r(J));
        %---1
        fprintf(fid,'    texture{TdG_LightSand}');
%         fprintf(fid,'    texture{ crackle  scale 1.5 turbulence 0.1 \n');
%         fprintf(fid,'        texture_map {[0.00 Texture_W] [0.08 Texture_W] [0.08 Texture_S] [1.00 Texture_S]} \n');
%         fprintf(fid,'            scale 0.2 \n');
%         fprintf(fid,'        }\n');
        fprintf(fid,'    matrix <   %12.4e,  %12.4e,  %12.4e, \n',matrix(J,1,:));
        fprintf(fid,'               %12.4e,  %12.4e,  %12.4e, \n',matrix(J,2,:));
        fprintf(fid,'               %12.4e,  %12.4e,  %12.4e, \n',matrix(J,3,:));
        fprintf(fid,'               %12.4e,  %12.4e,  %12.4e> \n',point(J,:));
        fprintf(fid,'        }\n');
    end
    %
    fprintf(fid,'mesh {\n');
    for J=1:1:length(bondedWallMesh)
        fprintf(fid,'triangle{\n');
        fprintf(fid,'<%25.16e,%25.16e,%25.16e>,\n',bondedWallPoint(bondedWallMesh(J,1)+1,:));
        fprintf(fid,'<%25.16e,%25.16e,%25.16e>,\n',bondedWallPoint(bondedWallMesh(J,2)+1,:));
        fprintf(fid,'<%25.16e,%25.16e,%25.16e>',bondedWallPoint(bondedWallMesh(J,3)+1,:));
        fprintf(fid,'}\n');
    end
    %
    fprintf(fid,'texture {T_Wood1}\n');
    %
    fprintf(fid,'scale <1,1,1>\n');
    fprintf(fid,'matrix <%25.16e,%25.16e,%25.16e, \n', bondedWallMatrix(1,1),bondedWallMatrix(1,2),bondedWallMatrix(1,3));
    fprintf(fid,'        %25.16e,%25.16e,%25.16e, \n', bondedWallMatrix(2,1),bondedWallMatrix(2,2),bondedWallMatrix(2,3));
    fprintf(fid,'        %25.16e,%25.16e,%25.16e, \n', bondedWallMatrix(3,1),bondedWallMatrix(3,2),bondedWallMatrix(3,3));
    %
    fprintf(fid,'        %25.16e,%25.16e,%25.16e> \n', [bondedWallX(I-999,1),bondedWallX(I-999,2),bondedWallX(I-999,3)]);  
    fprintf(fid,'}\n');
    %
    fprintf(fid,'cylinder{<2.5,0.0,0.0>,<-2.5,0.0,0.0>,1.0 texture {T_Wood1}}\n');
end


