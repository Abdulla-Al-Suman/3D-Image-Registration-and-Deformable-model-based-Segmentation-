clc;
clear all
close all

load('Patient 25 original MRI.mat')
I0 = imresize(CT,0.25);
[M,N,OO] = size(I0);
load('Patient 41 original MRI.mat')
R0 = imresize(CT,0.25);

I0 = I0-min(min(min(I0)));
I0 = I0.*255./max(max(max(I0)));
R0 = R0-min(min(min(R0)));
R0 = R0.*255./max(max(max(R0)));

p=13;

h1=fspecial('gauss',[p p],0.1);

h2=fspecial('gauss',[p p],1);

h3=fspecial('gauss',[p p],2);


for i=1:OO
    t1data(:,:,i)=filter2(h1,squeeze(R0(:,:,i)));
    
    t2data(:,:,i)=filter2(h1,squeeze(I0(:,:,i)));
    
end

% figure(3), imagesc(sum(t1data,3))

for i=1:OO
    t1data1(:,:,i)=filter2(h2,squeeze(R0(:,:,i)));
    t2data1(:,:,i)=filter2(h2,squeeze(I0(:,:,i)));
%     t1data1(:,:,i)=filter2(h2,squeeze(t1data(:,:,i)));
%     t2data1(:,:,i)=filter2(h2,squeeze(t2data(:,:,i)));
%     
end

for i=1:OO
    t1data2(:,:,i)=filter2(h3,squeeze(R0(:,:,i)));   
    t2data2(:,:,i)=filter2(h3,squeeze(I0(:,:,i)));
%     t1data2(:,:,i)=filter2(h3,squeeze(t1data1(:,:,i)));   
%     t2data2(:,:,i)=filter2(h3,squeeze(t2data1(:,:,i)));
    
end

T11{1}=t1data2;
T11{2}=t1data1;
T11{3}=t1data;

T22{1}=t2data2;
T22{2}=t2data1;
T22{3}=t2data;

save registration-input-p41T11P25T22 T11 T22

