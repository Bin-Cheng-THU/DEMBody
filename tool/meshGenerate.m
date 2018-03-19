clc;clear;
dx = 2.5; %2.5Rmax
x = 75; %Boundary size
y = 75; %Boundary size
z = 0;
% 为了防止颗粒一直位于网格边界上
% 出现数值不稳定的现象（即两个网格）
% 来回跑，因此将其设于网格内部
% 因此需要长方体偏心网格
% 上述说明错误！！！！！！

num1 = x/dx*2;
num2 = num1*num1;

casenum(1) = 1;
casenum(2) = -num1-1;
casenum(3) = -num1;
casenum(4) = -num1+1;
casenum(5) = -num2;
casenum(6) = -num2+1;
casenum(7) = -num2+num1-1;
casenum(8) = -num2+num1;
casenum(9) = -num2+num1+1;
casenum(10) = -num2-1;
casenum(11) = -num2-num1+1;
casenum(12) = -num2-num1;
casenum(13) = -num2-num1-1;
casenum = casenum';
disp(round(casenum(1)))
disp(round(casenum(2)))
disp(round(casenum(3)))
disp(round(casenum(4)))
disp(round(casenum(5)))
disp(round(casenum(6)))
disp(round(casenum(7)))
disp(round(casenum(8)))
disp(round(casenum(9)))
disp(round(casenum(10)))
disp(round(casenum(11)))
disp(round(casenum(12)))
disp(round(casenum(13)))

