%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Demo texture for POV-Ray
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all;
format long;
%% Head file
%---camera settings
bgclr = [0.05, 0.05, 0.05];
angle = 4.5;
location = [100.0, 0.0, 50.0];
skyvec = [0,0,1];
focus = [20.0, 0.0, 40.0];
up = [0.0, 0.9, 0.0];
right = [1.6, 0.0, 0.0];
lightsrc = [15.0, 15.0, 300.0];

fid = fopen(Demo,'wt');
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
%% Texture for Asteroid Regolith
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
fprintf(fid,'#declare asteroid_regolith_1= \n');
fprintf(fid,'  texture {pigment { White_Marble} finish {ambient 0.0 diffuse 0.7 brilliance 2.0 specular 0.05 roughness 0.1} normal { agate 0.13 scale 0.08 }} \n');
%% Texture for Asteroid Regolith
fprintf(fid,'#declare asteroid_regolith_2= \n');
fprintf(fid,'  texture   { pigment {  crackle form <1.0,0,0>  color_map { [0 rgb <0.15,0.15,0.15>] [1 rgb <0.3,0.3,0.3>] } cubic_wave  scale 0.1 \n');
fprintf(fid,'                         warp { turbulence .5  octaves 3  omega 1.0 lambda .7 } scale 0.5 } \n');
fprintf(fid,'              finish {ambient 0.0 diffuse 1.0 brilliance 2.0 specular 0.05 roughness 0.1 } \n');
fprintf(fid,'              normal { agate 0.13 scale 0.08 }\n');
fprintf(fid,'} \n');
%% Texture for Robot
fprintf(fid,'#declare PaintColor = color White; \n');
fprintf(fid,'#declare PaintBright = pigment{PaintColor} \n');
fprintf(fid,'#declare PaintDark = pigment{PaintColor / 2} \n');
fprintf(fid,'#declare CarPaint =  \n');
fprintf(fid,'  texture{  pigment{   aoi  pigment_map{ [0.5 PaintDark] [1.0 PaintBright]  } } \n');
fprintf(fid,'    normal {bozo 0.05 scale 0.1} \n');
fprintf(fid,'    finish{ diffuse 0.65 brilliance 0.6  reflection{ rgb <.05, .05, .05>, rgb<.2,.2,.2> fresnel on } \n');
fprintf(fid,'    } \n');
fprintf(fid,'  } \n');
%% Texture for Rubber
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
%% 
fprintf(fid,'#declare Sand= \n');
fprintf(fid,'  texture{  pigment{ color rgb <.518, .339, .138> } \n');
fprintf(fid,'    normal{ bumps 5 scale 0.05 } \n');
fprintf(fid,'    finish{ specular .3 roughness .8 } \n');
fprintf(fid,'  } \n');
fprintf(fid,'  texture{ pigment{ wrinkles scale 0.05 color_map{ \n');
fprintf(fid,'	[0.0 color rgbt <1, .847, .644, 0>] \n');
fprintf(fid,'	[0.2 color rgbt <.658, .456, .270, 1>] \n');
fprintf(fid,'	[0.4 color rgbt <.270, .191, .067, .25>] \n');
fprintf(fid,'	[0.6 color rgbt <.947, .723, .468, 0>] \n');
fprintf(fid,'	[0.8 color rgbt <.356, .250, .047, 1>] \n');
fprintf(fid,'	[1.0 color rgbt <.171, .136, .1, 1>] \n');
fprintf(fid,'      } \n');
fprintf(fid,'    } \n');
fprintf(fid,'  } \n');