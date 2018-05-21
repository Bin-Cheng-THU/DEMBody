%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           convert VTK files to POV-Ray files
%           input: mesh file; point file; bondedWall file; Foot file
%           output: POV-Ray files
%
%           description: Hopping Robot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear
format long;
%% Parameters
%---folder path of particle datas
folderParticle = 'D:\0-5--hopping robot\Data\';
folderPOVRay = 'D:\0-5--hopping robot\Data\';
%---number of files
Num = 100;
Step = 1;
TagQ = 0; % whether use Quanternion of particles
%---color setting
min = 0.0;
max = 1.0;
TagColor = 0; % whether use Color Map of particles
%---camera settings
bgclr = [0.05, 0.05, 0.05];
angle = 4.5;
location = [100.0, 0.0, 50.0];
skyvec = [0,0,1];
focus = [20.0, 0.0, 40.0];
up = [0.0, 0.9, 0.0];
right = [1.6, 0.0, 0.0];
lightsrc = [15.0, 15.0, 300.0];
%% Read
filename = strcat(folderParticle,'BondedWalls.txt');
tmp = importdata(filename);
bondedWallX = tmp(1:Step:end,[2 3 4]);
bondedWallQ = tmp(1:Step:end,[11 12 13 14]);
clear tmp;

filename = strcat(folderParticle,'BondedWallMesh.txt');
tmp = importdata(filename);
bondedWallMesh = tmp(:,[2 3 4]);
clear tmp;

filename = strcat(folderParticle,'BondedWallPoint.txt');
tmp = importdata(filename);
bondedWallPoint = tmp;
clear tmp;

filename = strcat(folderParticle,'FootMesh.txt');
tmp = importdata(filename);
FootMesh = tmp(:,[2 3 4]);
clear tmp;

filename = strcat(folderParticle,'FootPoint.txt');
tmp = importdata(filename);
FootPoint = tmp;
clear tmp;
%% Pov-Ray file generating
for I = 1000:(1000+Num-1)
    disp(I)
    %% read particle data
    filename = strcat(folderParticle,num2str(I));
    filename = strcat(filename,'X.csv');
    
    tmp = importdata(filename);
    tmp = tmp.data;
    point = tmp(:,1:3);
    r = tmp(:,10);
    if (TagQ)
        q = tmp(:,19:22);
    end
    if (TagColor)
        velocity = tmp(:,4:6);
        vMag = sqrt(velocity(:,1).^2 + velocity(:,2).^2 + velocity(:,3).^2);
    end
    clear tmp;
    
    %% attitude matrix of particles
    if (TagQ)
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
    end
    
    %% attitude matrix of bonded walls
    bondedWallMatrix = quat2dcm([bondedWallQ(I-999,4) bondedWallQ(I-999,1) bondedWallQ(I-999,2) bondedWallQ(I-999,3)]);
   
    %% write files
    filename = strcat(folderPOVRay,'Scene');
    filename = strcat(filename,num2str(I));
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
    fprintf(fid,'light_source {<%12.4e,%12.4e,%12.4e> color 2.5}\n',lightsrc);
    fprintf(fid,'light_source {<%12.4e,%12.4e,%12.4e> color 0.5}\n',location);
    fprintf(fid,'background { color rgb <%12.4e, %12.4e, %12.4e> }\n', bgclr);
    fprintf(fid,'camera {location <%12.4e,%12.4e,%12.4e> \n', location);
    fprintf(fid,'        sky   <%12.4e,%12.4e,%12.4e> \n', skyvec);
    fprintf(fid,'        look_at <%12.4e,%12.4e,%12.4e> \n', focus);
    fprintf(fid,'        up <%12.4e,%12.4e,%12.4e> \n', up);
    fprintf(fid,'        right <%12.4e,%12.4e,%12.4e> \n', right);
    fprintf(fid,'        }\n'); 
    %
    fprintf(fid,'#declare White_Marble_Map = \n');
    fprintf(fid,'color_map { \n');
    fprintf(fid,'    [0.0 rgb <0.5, 0.5, 0.5>] \n');
    fprintf(fid,'    [0.8 rgb <0.3, 0.3, 0.3>] \n');
    fprintf(fid,'    [1.0 rgb <0.1, 0.1, 0.1>] \n');
    fprintf(fid,'} \n');
    fprintf(fid,'#declare White_Marble = \n');
    fprintf(fid,'pigment { \n');
    fprintf(fid,'    marble    turbulence 1.0    color_map { White_Marble_Map } \n');
    fprintf(fid,'} \n');
    fprintf(fid,'#declare asteroid_regolith= \n');
    fprintf(fid,'  texture {pigment { White_Marble} finish {ambient 0.0 diffuse 0.7 brilliance 2.0 specular 0.05 roughness 0.1} normal { agate 0.13 scale 0.08 }} \n');
    %
    fprintf(fid,'#declare Rubber = \n');
    fprintf(fid,'  texture{ pigment{ aoi color_map{ \n');
    fprintf(fid,'	[0.00 rgb <.0075, .0175, .0025>] \n');
    fprintf(fid,'	[0.55 rgb <.020, .022, .024>] \n');
    fprintf(fid,'	[0.65 rgb <.004, .004, .004>] \n');
    fprintf(fid,'	[0.85 rgb <.006, .002, .001>] \n');
    fprintf(fid,'	[1.00 rgb <.007, .004, .001>]} \n');
    fprintf(fid,'    poly_wave 1.25  scale 0.1  } \n');
    fprintf(fid,'    normal {bozo 0.2 scale 0.25} \n');
    fprintf(fid,'    finish{ specular .015  roughness .075  brilliance 0.275 \n');
    fprintf(fid,'    } \n');
    fprintf(fid,'} \n');
    
    num = 0;
    for J = 1:length(r)
        %-------------
        Flag = 1;
        if (and(and(abs(point(J,1))<50,abs(point(J,2))<50),abs(point(J,3)-15)<7))
            Flag = 0;
        end
        
        if (Flag==1)
            num = num + 1;
            if (TagQ)
                fprintf(fid,'sphere ');
                fprintf(fid,'{<%25.16e,%25.16e,%25.16e>, %25.16e \n',[0,0,0],1);
                fprintf(fid,'  texture{asteroid_regolith} \n');
                fprintf(fid,'  matrix <   %12.4e,  %12.4e,  %12.4e, \n',matrix(J,1,:));
                fprintf(fid,'             %12.4e,  %12.4e,  %12.4e, \n',matrix(J,2,:));
                fprintf(fid,'             %12.4e,  %12.4e,  %12.4e, \n',matrix(J,3,:));
                fprintf(fid,'             %12.4e,  %12.4e,  %12.4e> \n',point(J,:));
                fprintf(fid,'  scale %12.4e \n',r(J));
                fprintf(fid,'} \n');
            else          
                fprintf(fid,'sphere ');
                fprintf(fid,'{<%25.16e,%25.16e,%25.16e>, %25.16e \n',point(J,:),r(J));
                fprintf(fid,'    texture{asteroid_regolith} \n');
                fprintf(fid,'        }\n'); 
            end
        end
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
    fprintf(fid,'texture {pigment {color rgb<0.625000, 0.425000, 0.025000> } finish {reflection 0.08}} \n');
    fprintf(fid,'scale <1,1,1>\n');
    fprintf(fid,'matrix <%25.16e,%25.16e,%25.16e, \n', bondedWallMatrix(1,1),bondedWallMatrix(1,2),bondedWallMatrix(1,3));
    fprintf(fid,'        %25.16e,%25.16e,%25.16e, \n', bondedWallMatrix(2,1),bondedWallMatrix(2,2),bondedWallMatrix(2,3));
    fprintf(fid,'        %25.16e,%25.16e,%25.16e, \n', bondedWallMatrix(3,1),bondedWallMatrix(3,2),bondedWallMatrix(3,3));
    fprintf(fid,'        %25.16e,%25.16e,%25.16e> \n', [bondedWallX(I-999,1),bondedWallX(I-999,2),bondedWallX(I-999,3)]);  
    fprintf(fid,'}\n');    
    %
    fprintf(fid,'mesh {\n');
    for J=1:1:length(FootMesh)
        fprintf(fid,'triangle{\n');
        fprintf(fid,'<%25.16e,%25.16e,%25.16e>,\n',FootPoint(FootMesh(J,1)+1,:));
        fprintf(fid,'<%25.16e,%25.16e,%25.16e>,\n',FootPoint(FootMesh(J,2)+1,:));
        fprintf(fid,'<%25.16e,%25.16e,%25.16e>',FootPoint(FootMesh(J,3)+1,:));
        fprintf(fid,'}\n');
    end
    fprintf(fid,'texture  {Rubber} \n');
    fprintf(fid,'scale <1,1,1>\n');
    fprintf(fid,'matrix <%25.16e,%25.16e,%25.16e, \n', bondedWallMatrix(1,1),bondedWallMatrix(1,2),bondedWallMatrix(1,3));
    fprintf(fid,'        %25.16e,%25.16e,%25.16e, \n', bondedWallMatrix(2,1),bondedWallMatrix(2,2),bondedWallMatrix(2,3));
    fprintf(fid,'        %25.16e,%25.16e,%25.16e, \n', bondedWallMatrix(3,1),bondedWallMatrix(3,2),bondedWallMatrix(3,3));
    fprintf(fid,'        %25.16e,%25.16e,%25.16e> \n', [bondedWallX(I-999,1),bondedWallX(I-999,2),bondedWallX(I-999,3)]);  
    fprintf(fid,'}\n');
end
fclose('all');


