%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           combine Movie
%           input: several movies
%           output: Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
format long;
%% parameters
moviePath = 'C:\Users\chengbin\Desktop\';

%%
filename = strcat(moviePath,'out1.avi');
aviobj1 = VideoReader(filename);

filename = strcat(moviePath,'out2.avi');
aviobj2 = VideoReader(filename);

%%
fps = aviobj1.FrameRate;
% N = aviobj1.NumberOfFrames + aviobj2.NumberOfFrames;

aviobj = VideoWriter('out.avi'); %初始化一个avi文件
fps = aviobj.FrameRate;
open(aviobj);

while hasFrame(aviobj1)
    frames = readFrame(aviobj1); %avi需要彩色数据
    writeVideo(aviobj,frames); %一帧一帧的写入avi
end

while hasFrame(aviobj2)
    frames = readFrame(aviobj2); %avi需要彩色数据
    writeVideo(aviobj,frames); %一帧一帧的写入avi
end
close(aviobj); %将缓存数据写入avi
