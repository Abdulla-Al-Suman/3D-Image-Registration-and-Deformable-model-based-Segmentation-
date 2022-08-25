function run_range_test
clc
clear all
close all
global PDF num_iter filter_params sm Im V VV MFEF x y
tic
PDF = 15;% PDF = parameter dividing factor
s = 3; %s= parameter selecting number
sm = 0; %sm= similarity measurement, 0= SCV,1=MI
[x,y] = meshgrid((0:128-1),(0:128-1));
% number of iterations at each level
num_iter = [10 5 5];
% size and variance of laplacian-of-gaussian filters for each level
filter_params = [10 0.9;5 0.4;3 0.02];
load('elastic32VP41VVP25.mat')

m=[0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000];
%m=[0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000];
Im = m; % Im = Initial m for other slices
MFEF = zeros(45*4,18);


%for i = 1:1
    
[R0,I] = test_range(m,s);

V=I;VV=R0;MFEFL = MFEF;

save  elasticLVP41VVP25 V VV MFEFL
toc

%end

