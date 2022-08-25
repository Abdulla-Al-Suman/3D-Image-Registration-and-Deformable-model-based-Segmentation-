
%function run_range_test

clc;
close all
clear all;
global sm num_iter m x y z
sm = 0; %sm= similarity measurement, 0= SCV,1=MI
load('registration-input-p41T11P25T22.mat')
[x,y,z] = meshgrid([0:256-1]-256/2,[0:256-1]-256/2,[0:45-1]-45/2);
% number of iterations at each level
num_iter = [20 10 10];

m=[1.000    0.0000  0.0000  0.000   0.0000  1.000   0.0000  0.0000  0.0000  0.0000  1.0000  0.0000];
%m=[0.9745    0.0223   -0.1968   -4.0309  0.0033  0.8819  -0.2954  11.9348  0.0080  -0.0196  1.1436  -1.2354];
%m=[0.9113 0.0728 0.0131 7.1014 -0.0731 0.9112 0.0179 9.9215 -0.0117 -0.0189 0.9140 8.7884];
%m=[1.0896    -0.0869  -0.0134  -6.7514   0.0867  1.0899   -0.0221  -11.2346  0.0152  0.0209  1.0929  -9.9326];
%m=[0.9113    0.0728  0.0131  7.1014   -0.0731  0.9112   0.0179  9.9215  -0.0117  -0.0189  0.9140  8.7884];

% mex interp3_mex.cpp
% mex mul_hes_mex.cpp                                                                  

tic 
[R0,I] = test_range(T11,T22);
figure(5), imagesc([sum(R0,3)   sum(I,3)]), colormap(gray)    
toc
m
V=I;VV=R0;
save affine-NI40iVP41VVP25 V VV m
