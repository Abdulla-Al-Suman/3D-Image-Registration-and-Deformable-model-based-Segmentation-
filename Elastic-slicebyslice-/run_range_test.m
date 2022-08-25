function run_range_test
clc
clear all
close all

global PDF num_iter filter_params sm Im V VV MFEF x y
PDF = 100;% PDF = parameter dividing factor
s = 4; %s= parameter selecting number
sm = 0; %sm= similarity measurement, 0= SCV,1=MI 
load('elastic18VP41VVP25.mat') % Contains V VV
[M,N,OO] = size(V);
[x,y] = meshgrid((0:256-1),(0:256-1));
%[x,y] = meshgrid([0:256-1]-256/2,[0:256-1]-256/2);
%x = 3*ones(256,256); y = 3*ones(256,256);
% number of iterations at each level
num_iter = [20 20 20];

% size and variance of laplacian-of-gaussian filters for each level
filter_params = [30 2.02;15 1.5;10 0.08];

MFEF = zeros(OO,32);% MFEF = m for each frame

m=[0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000];
%m=[0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000];
%m=[0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000];
Im = m; % Im = Initial m for other slices

%for i = 1:1  % this loop for using many parameter, using rand_start_posns_uniform, for registration program
tic    
[R0,I] = test_range(m,s);
toc
V=I;VV=R0; MFEF32 = MFEF;
save  elastic32VP41VVP25 V VV MFEF32

%end

