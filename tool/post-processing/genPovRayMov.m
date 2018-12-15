%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           convert POV-Ray files to Movie
%           input: POV-Ray files
%           output: Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
format long;
%% parameters
fps = 10;
N = 193;
folderPath = 'C:\Users\chengbin\Desktop\genWallPovRay\genWallPovRay\Data\';
aviobj = VideoWriter('out-yorp-x.avi'); %初始化一个avi文件
aviobj.FrameRate = fps;
open(aviobj);

for ii=1000:(1000+N-1)
    filename = strcat(folderPath,'Scene');
    filename = strcat(filename,num2str(ii));
    filename = strcat(filename,'.png');
    
    frames=imread(filename); %avi需要彩色数据
    writeVideo(aviobj,frames); %一帧一帧的写入avi
end
close(aviobj); %将缓存数据写入avi