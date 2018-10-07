clc;clear;
format long;
%%
check = zeros(202,3);
for ii = 1:length(check)
    ii
    filename = num2str((ii-1)*0.0001,'%10.5f');
    filename = [filename,'Force.txt'];
    file1 = ['Lattice\',filename];
    file2 = ['All\',filename];
    data1 = load(file1);
    data2 = load(file2);
%     data1 = sortrows(sortrows(data1,1),2);
%     data2 = sortrows(sortrows(data2,1),2);
    error = abs(data1-data2);
    check(ii,:) = sum(error);
end