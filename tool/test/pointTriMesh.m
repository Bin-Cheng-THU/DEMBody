% Check which part of a triangle a points is.
clc;clear;
format long;

%% Input parameters
global pointA pointB pointC;
pointA = [-1   0   0];
vectorN = [0   0   1];
point = [2.08937093999600 8.74381155525819 4.62698833028458];
RVb = point-pointA-dot(point-pointA,vectorN)*vectorN;
%% 
AB = [2  1  0];%pointB - pointA;
AC = [0  2  0];%pointC - pointA;
pointB = pointA+AB;
pointC = pointA+AC;
AP = RVb;

% u = (-dot(AP,AB)*dot(AB,AC)+dot(AP,AC)*dot(AB,AB))/(-dot(AC,AB)*dot(AB,AC)+dot(AC,AC)*dot(AB,AB));
% v = (dot(AP,AB)*dot(AC,AC)-dot(AP,AC)*dot(AC,AB))/(dot(AB,AB)*dot(AC,AC)-dot(AB,AC)*dot(AC,AB));
% 
% un = (dot(AP,AB)-dot(AP,AC)-dot(AB,AB)+dot(AB,AC))/(2*dot(AC,AB)-dot(AC,AC)-dot(AB,AB));
% vn = 1-un;
% pointn = un*AC+vn*AB+pointA;
% 
Tag = Check(pointA,pointB,pointC,RVb+pointA);
disp(Tag);

%% Plot
figure (1)
hold on;
axis equal;
plot3([pointA(1),pointB(1),pointC(1),pointA(1)],[pointA(2),pointB(2),pointC(2),pointA(2)],[pointA(3),pointB(3),pointC(3),pointA(3)]);
for ii = 1:5000
    u = rand()*4-2;
    v = rand()*4-2;
    point = u*AC + v*AB + pointA;
    Tag = Check(pointA,pointB,pointC,point);
    switch Tag
        case {1}
            color = 'k';
        case {2}
            color = 'b';
        case {3}
            color = 'y';
        case {4}
            color = 'g';
        case {5}
            color = 'r';
        case {6}
            color = 'm';
        case {7}
            color = 'c';
        case {8}
            color = 'k';
    end
    
    switch Tag
        case {1,2,3,4,5,6,7}
            size = 2;
        case {8}
            size = 0.5;
    end
    plot3(point(1),point(2),point(3),'Marker','o','MarkerSize',size,'Color',color);
end

