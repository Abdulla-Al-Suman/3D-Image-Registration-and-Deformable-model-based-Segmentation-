% Registered data visualization
clc
clear all 
close all
load('Patient 25 muscle data.mat')
load('elastic18VP41VVP25.mat')
x1=V;x2=VV;
[m,n,z] = size(x1);

for i=10:30
% Retyriving x y data for left side muscles  
xl3 = muscle_data{3}{1}{i}.x;yl3 = muscle_data{3}{1}{i}.y;
xl11 = muscle_data{11}{1}{i}.x;yl11 = muscle_data{11}{1}{i}.y;
xl16 = muscle_data{16}{1}{i}.x;yl16 = muscle_data{16}{1}{i}.y;
xl22 = muscle_data{22}{1}{i}.x;yl22 = muscle_data{22}{1}{i}.y;
xrl3 = xl3+256;xrl11 = xl11+256;xrl16 = xl16+256; xrl22 = xl22+256;

% Retyriving x y data for right side muscles

xr3 = muscle_data{3}{2}{i}.x;yr3 = muscle_data{3}{2}{i}.y;
xr11 = muscle_data{11}{2}{i}.x;yr11 = muscle_data{11}{2}{i}.y;
xr16 = muscle_data{16}{2}{i}.x;yr16 = muscle_data{16}{2}{i}.y;
xr22 = muscle_data{22}{2}{i}.x;yr22 = muscle_data{22}{2}{i}.y;
xrr3 = xr3+256;xrr11 = xr11+256; xrr16 = xr16+256;xrr22 = xr22+256;
%%%%%%%%%

imagesc([x2(:,:,i) x1(:,:,i)]), colormap(gray(256));
hold on
plot(xl3,yl3,'g',xl11,yl11,'g',xl16,yl16,'r',xl22,yl22,'c',xr3,yr3,'b',xr11,yr11,'b',xr16,yr16,'m',xr22,yr22,'b')
plot(xrl3,yl3,'g',xrl11,yl11,'g',xrl16,yl16,'r',xrl22,yl22,'b',xrr3,yr3,'c',xrr11,yr11,'c',xrr16,yr16,'m',xrr22,yr22,'c')
pause
    
end


