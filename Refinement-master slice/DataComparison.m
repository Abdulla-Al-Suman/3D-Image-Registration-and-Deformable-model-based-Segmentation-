clc
clear all 
close all
% load('registeredimage_whiplash_t1t2.mat')
% load('registeredimage_whiplash.mat')
%load('elastic_output8.mat')
load('p2 muscle data.mat')
load('FFR81832p-808060m104050m-1.mat')
 x1=V;
 x2=VV;
[m,n,z] = size(x1);

% [m,n,z] = size(x1);% load('scv_m_data_results.mat')
% x1=I;
% x2=R0;
% [m,n,z] = size(x1);









% load('Patient2_11.mat')
% x1=V;
% [m,n,z] = size(x1);
% load('P2_16.mat')
% x2=V;
% for i=1:256
%     
%     ct1=squeeze(x1(:,i,:));
%     ct2=squeeze(x2(:,i,:));
%     imagesc([ct1 ct2]), colormap(gray(256));
%     
%     pause
%     
% end
for i=1:z
  % Retyriving x y data for left side muscles  
xl1 = muscle_data{1}{1}{i}.x;
yl1 = muscle_data{1}{1}{i}.y;
xl2 = muscle_data{2}{1}{i}.x;
yl2 = muscle_data{2}{1}{i}.y;
% xl3 = muscle_data{3}{1}{z}.x;
% yl3 = muscle_data{3}{1}{z}.y;
% xl4 = muscle_data{4}{1}{z}.x;
% yl4 = muscle_data{4}{1}{z}.y;
xl5 = muscle_data{5}{1}{i}.x; % centre
yl5 = muscle_data{5}{1}{i}.y;
% xl7 = muscle_data{7}{1}{z}.x;
% yl7 = muscle_data{7}{1}{z}.y;
xrl1 = xl1+256;xrl5 = xl5+256; xrl2 = xl2+256;%xrl3 = xl3+256; xrl4 = xl4+256; xrl7 = xl7+256;


% Retyriving x y data for right side muscles

xr1 = muscle_data{1}{2}{i}.x;
yr1 = muscle_data{1}{2}{i}.y;
xr2 = muscle_data{2}{2}{i}.x;
yr2 = muscle_data{2}{2}{i}.y;
% xr3 = muscle_data{3}{2}{z}.x;
% yr3 = muscle_data{3}{2}{z}.y;
% xr4 = muscle_data{4}{2}{z}.x;
% yr4 = muscle_data{4}{2}{z}.y;
% xr5 = muscle_data{5}{2}{z}.x;
% yr5 = muscle_data{5}{2}{z}.y;
xrr1 = xr1+256; xrr2 = xr2+256;%xrr3 = xr3+256; xrr4 = xr4+256;xrr5 = xr5+256;
    %%%%%%%%%
    ct1=x1(:,:,i);
    ct2=x2(:,:,i);
%     imagesc([ct1 ct2]), colormap(gray(256));
imagesc([ct2 ct1]), colormap(gray(256));
hold on
plot(xl1,yl1,'g',xl5,yl5,'r',xr1,yr1,'b',xl2,yl2,'c',xr2,yr2,'m')
plot(xrl1,yl1,'g',xrl5,yl5,'r',xrr1,yr1,'b',xrl2,yl2,'c',xrr2,yr2,'m')

 
% imagesc(ct2); %colormap(gray(25);
%     
    pause
    
end


