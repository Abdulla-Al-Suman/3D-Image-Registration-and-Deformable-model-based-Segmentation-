% Original data visualization
clc
clear all 
close all
load('Patient 41 original MRI.mat');
x1= imresize(CT,0.25);[m,n,z] = size(x1);
load('Patient 25 original MRI.mat');x2= imresize(CT,0.25);
for i=1:z  
imagesc([x1(:,:,i) x2(:,:,i)]), colormap(gray(256));
pause   
end


