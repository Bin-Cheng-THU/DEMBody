clc;clear;
format long;
%%
check = zeros(65,3);
for ii = 1:length(check)
    ii
    filename = num2str((ii-1)*0.0002,'%10.5f');
    filename = [filename,'Force.txt'];
    file1 = ['1\',filename];
    file2 = ['3\',filename];
    data1 = load(file1);
    data2 = load(file2);
%     data1 = sortrows(sortrows(data1,1),2);
%     data2 = sortrows(sortrows(data2,1),2);
    error = abs(data1-data2);
    check(ii,:) = sum(error);
end