clc
clear all 
close all
% I(X,Y,Z) R(x,y,z)
I1 = [128;165;7]; I2 = [131;112;7]; I3 = [208;124;30];I4 = [94;103;30];
R1 = [134;151;11]; R2 = [132;105;11];R3 = [211;108;34];R4 = [103;99;34];
S = ((norm(I2-I1)/norm(R2-R1))+(norm(I3-I1)/norm(R3-R1))+(norm(I4-I1)/norm(R4-R1))+(norm(I3-I2)/norm(R3-R2))+(norm(I4-I2)/norm(R4-R2))+(norm(I4-I3)/norm(R4-R3)))/6;
%R1 = S*R1;R2 = S*R2;R3 = S*R3;R4 = S*R4;
I1 = I1/S;I2 = I2/S;I3 = I3/S;I4 = I4/S;
pprime = (I1+I2+I3+I4)/4;p = (R1+R2+R3+R4)/4;
q1prime = I1-pprime;q2prime = I2-pprime;q3prime = I3-pprime;q4prime = I4-pprime;
q1 = R1-p;q2 = R2-p;q3 = R3-p;q4 = R4-p; 
H = (q1*q1prime')+(q2*q2prime')+(q3*q3prime')+(q4*q4prime');
[U,SS,V] = svd(H);
X = V*U';
Deterrminant = det(X)
R = S*X
T = pprime-(X*p)





