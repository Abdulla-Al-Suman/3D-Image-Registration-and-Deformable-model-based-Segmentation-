clc
clear all 
close all

load('TMD-m22sLf12.mat')   % contains xnnnn, ynnnn for slice 25
load('P1.mat')
I = zeros(256,256);
I(:,:) = V(:,:,SNIRI);
load('p2 muscle data.mat')
MX = muscle_data{MN}{SNP}{SNIRI}.x; % manual data for comparison; 1=left
MY = muscle_data{MN}{SNP}{SNIRI}.y;
% X = muscle_data{5}{1}{3}.x;  % 5th muscle 1= left side, 1=slice x data
% Y = muscle_data{5}{1}{3}.y;
X = xnnnn;  % 5th muscle 1= left side, 1=slice x data
Y = ynnnn;
% BW = roipoly(I, X, Y);
% stats = regionprops(BW,'Centroid');
% a=stats.Centroid(1)'
% CNT =[stats.Centroid(1),stats.Centroid(2)]
% %CNT= BW.Centroid

%%%%%%%%%%%%%%%%%%
% x = 1:length(X);
% v = X;
% xq = 1:0.5:length(X);
% X=[];
% X = interp1(x,v,xq,'spline');                       %RESAMPLING
% 
% x = 1:length(Y);
% v = Y;
% xq = 1:0.5:length(Y);
% Y=[];
% Y = interp1(x,v,xq,'spline');
%%%%%%%%%%%%%%%%%%%%%%


G=fspecial('gauss',[40 40],3);
CI(:,:)=filter2(G,I);
[dCIdx,dCIdy] = gradient(CI);
dx = dCIdx.^2;
dy = dCIdy.^2;
GM = dx+dy;
[dCIdx1,dCIdy1] = gradient(GM);
GM1 = dCIdx1+dCIdy1;
%GM = (-13)*dCIdx+dCIdy*(-13);
%GM = 10*dCIdy;

figure
imagesc(GM), colormap(gray(256));
hold on
plot(X,Y,'b')
hold on
plot(MX,MY,'r')
% figure
% imagesc(GM), colormap(gray(256));
% figure
% imagesc(I), colormap(gray(256));