% Grouping Deformation Algorithm
%3 force act; internal(local curvature), Gaussian potential, external
%dynamic distance, 1st search along local curvature then along reverse
%local curvature, then calculate acceleration , velocity, position change,
%If current poins is boundary then no change, If don't get boundary point
%within search then acts two force(internal+external). filter applied in
%four stages varing sigma. Remove overlap in each iteration for each muscle

clc
clear all 
close all
% Loading transfered labels
load('TMD-m3sLf12.mat')% contains xnnnn, ynnnn,xnnnnn,ynnnnn
X3 = xnnnnn; % X3=3 no muscle
Y3= ynnnnn; load('TMD-m16sLf12.mat');X16 = xnnnnn; Y16= ynnnnn; load('TMD-m19sLf12.mat');X19 = xnnnnn; Y19= ynnnnn;
load('TMD-m22sLf12.mat');X22 = xnnnnn; Y22= ynnnnn;

% Loading Nobel image volume and ground truth
load('P1.mat')
load('p2 muscle data.mat') %p2=P1, contains manual contour , can be used for comparison
MX3 = muscle_data{3}{SNP}{SNIRI}.x; % manual data for comparison; 1=left
MY3 = muscle_data{3}{SNP}{SNIRI}.y;
MX16 = muscle_data{16}{SNP}{SNIRI}.x; MY16 = muscle_data{16}{SNP}{SNIRI}.y;MX19 = muscle_data{19}{SNP}{SNIRI}.x; MY19 = muscle_data{19}{SNP}{SNIRI}.y;
MX22 = muscle_data{22}{SNP}{SNIRI}.x; MY22 = muscle_data{22}{SNP}{SNIRI}.y;
I = zeros(256,256);
I(:,:) = V(:,:,SNIRI);  % Nobel Image

figure(1)
imagesc(I), colormap(gray(256));
hold on
plot(X3,Y3,'b',MX3,MY3,'g',X16,Y16,'r',MX16,MY16,'g',X19,Y19,'c',MX19,MY19,'g',X22,Y22,'y',MX22,MY22,'g')

%%%%%%%%%%%%%%%%%%
x = 1:length(X3);v = X3;xq = 1:0.5:length(X3);X3=[];
X3 = interp1(x,v,xq,'spline');                       %RESAMPLING muscle3
x = 1:length(Y3);v = Y3;xq = 1:0.5:length(Y3);Y3=[];
Y3 = interp1(x,v,xq,'spline');
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
x = 1:length(X16);v = X16;xq = 1:0.5:length(X16);X16=[];
X16 = interp1(x,v,xq,'spline');                       %RESAMPLING muscle16
x = 1:length(Y16);v = Y16;xq = 1:0.5:length(Y16);Y16=[];
Y16 = interp1(x,v,xq,'spline');
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
x = 1:length(X19);v = X19;xq = 1:0.5:length(X19);X19=[];
X19 = interp1(x,v,xq,'spline');                       %RESAMPLING muscle19
x = 1:length(Y19);v = Y19;xq = 1:0.5:length(Y19);Y19=[];
Y19 = interp1(x,v,xq,'spline');
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
x = 1:length(X22);v = X22;xq = 1:0.5:length(X22);X22=[];
X22 = interp1(x,v,xq,'spline');                       %RESAMPLING muscle22
x = 1:length(Y22);v = Y22;xq = 1:0.5:length(Y22);Y22=[];
Y22 = interp1(x,v,xq,'spline');
%%%%%%%%%%%%%%%%%%%%%%
% saving initial positions
XXX3 = X3;
YYY3 = Y3; 
XXX16 = X16;XXX19 = X19;XXX22 = X22;YYY16 = Y16;YYY19 = Y19;YYY22 = Y22;

%Centroid Calculation muscle3
ROIM3 = roipoly(I, X3, Y3); %ROIM =ROI Musk
CENS3 = regionprops(ROIM3,'Centroid'); %CENS =Centroid Structure
CNTy3 =  CENS3.Centroid(2)-0.5;CNTy3 = -(CNTy3);
CNT3 =[CENS3.Centroid(1)-0.5,CNTy3]; %CNT=Centroid

%Centroid Calculation muscle16
ROIM16 = roipoly(I, X16, Y16); CENS16 = regionprops(ROIM16,'Centroid'); CNTy16 =  CENS16.Centroid(2)-0.5;
CNTy16 = -(CNTy16);CNT16 =[CENS16.Centroid(1)-0.5,CNTy16];
%Centroid Calculation muscle19
ROIM19 = roipoly(I, X19, Y19); CENS19 = regionprops(ROIM19,'Centroid'); CNTy19 =  CENS19.Centroid(2)-0.5;
CNTy19 = -(CNTy19);CNT19 =[CENS19.Centroid(1)-0.5,CNTy19];
%Centroid Calculation muscle22
ROIM22 = roipoly(I, X22, Y22); CENS22 = regionprops(ROIM22,'Centroid'); CNTy22 =  CENS22.Centroid(2)-0.5;
CNTy22 = -(CNTy22);CNT22 =[CENS22.Centroid(1)-0.5,CNTy22];
% ROI= or(ROIM3,ROIM16);
% ROI1= or(ROI,ROIM19);
% figure
% imagesc(ROI1), colormap(gray(256));

% Transfering contour into cartesian coordinate
X3=X3-0.5;Y3=Y3-0.5;Y3 = -(Y3);L3 = length(X3);X16=X16-0.5;Y16=Y16-0.5;Y16 = -(Y16);L16 = length(X16);
X19=X19-0.5;Y19=Y19-0.5;Y19 = -(Y19);L19 = length(X19);X22=X22-0.5;Y22=Y22-0.5;Y22 = -(Y22);L22 = length(X22);

figure(2)
imagesc(I), colormap(gray(256));
hold on
plot(X3+0.5,-(Y3)+0.5,'b',X16+0.5,-(Y16)+0.5,'r',X19+0.5,-(Y19)+0.5,'g',X22+0.5,-(Y22)+0.5,'y')

sigma=9;TGM3 =1000;TGM16 =100;TGM19 =1000;TGM22 =150; Dmax3 = 50;Dmax16 = 20;Dmax19 = 20;Dmax22 = 15;
win3 = 10;win16 = 10;win19 = 1000;win22 = 5; % internal force waiting
wex3 = 2;wex16 = 50;wex19 = 50;wex22 = 100;wd3 = 10;wd16 = 17.8;wd19 = 17.8;wd22 = 2; %im=20;%wdamp=1.6;
        %TSDmax = 0;
for isd=1:4 %isd=iteration starnard diviation;  filter varing loop
            %isd
   k1 = [0 0 0 0 0 0 -.5 1 -.5 0 0 0 0 0 0];
   win3 = win3/10000;win16 = win16/10000;win19 = win19/10000;win22 = win22/10000; wex3 = wex3/10000;wex16 = wex16/10000;wex19 = wex19/10000;wex22 = wex22/10000;
   wd3 = wd3/3;wd16 = wd16/3;wd19 = wd19/3;wd22 = wd22/3; wdamp = -.8;  
   Dmax3 = ceil(Dmax3/3);Dmax16 = ceil(Dmax16/3);Dmax19 = ceil(Dmax19/3);Dmax22 = ceil(Dmax22/3); % distance threshold
   TGM3 =3*(isd)*TGM3;TGM16 =3*(isd)*TGM16;TGM19 =3*(isd)*TGM19;TGM22 =3*(isd)*TGM22; % TGM= threshold gradient magnitude
   sigma=sigma/3; m1 = 1; % m1= mass for all vertices
   %im=im/2;
   bp = [0 0]; % bp=boundary point
   Nv = [0;0]; % Nv=normal vector             
   SVMF=0; % stop vertex moving flag
   SD = 0;

  % Calculating potential energy function/distribution
  G=fspecial('gauss',[40 40],sigma);CI(:,:)=filter2(G,I);[dCIdx,dCIdy] = gradient(CI); dx = dCIdx.^2; dy = dCIdy.^2; GM = dx+dy;
  GM1 = dCIdx+dCIdy;[dCIdx1,dCIdy1] = gradient(GM1);GM2 = dCIdx1+dCIdy1;
  % G1=fspecial('gauss',[40 40],sigma);
  % CI1(:,:)=filter2(G1,I);[dx1,dy1] = gradient(CI1); dx1 = dx1.^2; dy1 = dy1.^2; GM = dx1+dy1;
  %GM = CI;
  % figure
  % imagesc(GM), colormap(gray(256));
  Eim = -(GM2); % - means valley
  [dEimdx,dEimdy] = gradient(Eim);
  dEimdx = -1*dEimdx; % external force field
  dEimdy = -1*dEimdy; 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Starting Deformation process
  XX3 = zeros(1,L3);XX16 = zeros(1,L16);XX19 = zeros(1,L19);XX22 = zeros(1,L22); % saving deformed contour x values for all vertices
  YY3 = zeros(1,L3);YY16 = zeros(1,L16);YY19 = zeros(1,L19);YY22 = zeros(1,L22);
  VVV3 = zeros(L3,2);VVV16 = zeros(L16,2);VVV19 = zeros(L19,2);VVV22 = zeros(L22,2); % initial velocities       
  %SV = zeros(1,Dmax); % SV=search vector        

       for i=1:10 % deformation iteration
                    %i
                       NearestMR3 = roipoly(I, X22+0.5, -(Y22)+0.5); % NearestMR3 = Neatrest muscle region for m3 creating binary mask, 1 inside                     
                       % Checking m3 overlap
                           for overlapj3=1:L3  %overlapj=overlap checking j                              
                                    if NearestMR3(round(abs(Y3(overlapj3))+0.5),round(X3(overlapj3)+0.5))==1
                                      [Ny3,Nx3]= removeoverlap(overlapj3,X3,Y3,L3,CNT3,NearestMR3);
                                      X3(overlapj3) = Nx3; Y3(overlapj3) = Ny3;
                                    end
                           end
                       [X3,Y3] = deformation(TGM3,Dmax3,X3,Y3,XX3,YY3,VVV3,k1,dEimdx,dEimdy,GM,wd3,win3,wex3,wdamp,m1,L3,bp,SVMF,SD,Nv,NearestMR3); 
                     
                       N22MR16 = roipoly(I, X22+0.5, -(Y22)+0.5);N19MR16 = roipoly(I, X19+0.5, -(Y19)+0.5); %N22MR16= nearest 22 muscle region for 16
                       NearestMR16 = or(N22MR16,N19MR16);
                       figure(3)
                       imagesc(I), colormap(gray(256));
                       hold on
                       plot(X3+0.5,-(Y3)+0.5,'b',X16+0.5,-(Y16)+0.5,'r',X19+0.5,-(Y19)+0.5,'g',X22+0.5,-(Y22)+0.5,'y')
                          %pause
                          for overlapj16=1:L16  %overlapj=overlap checking j                              
                                    if NearestMR16(round(abs(Y16(overlapj16))+0.5),round((X16(overlapj16))+0.5))==1
                                      [Ny16,Nx16]= removeoverlap(overlapj16,X16,Y16,L16,CNT16,NearestMR16);
                                      X16(overlapj16) = Nx16; Y16(overlapj16) = Ny16;
                                    end
                          end
                       figure(4)
                       imagesc(I), colormap(gray(256));
                       hold on
                       plot(X3+0.5,-(Y3)+0.5,'b',X16+0.5,-(Y16)+0.5,'r',X19+0.5,-(Y19)+0.5,'g',X22+0.5,-(Y22)+0.5,'y')
                           %pause
                       [X16,Y16] = deformation(TGM16,Dmax16,X16,Y16,XX16,YY16,VVV16,k1,dEimdx,dEimdy,GM,wd16,win16,wex16,wdamp,m1,L16,bp,SVMF,SD,Nv,NearestMR16);

                       N22MR19 = roipoly(I, X22+0.5, -(Y22)+0.5);N16MR19 = roipoly(I, X16+0.5, -(Y16)+0.5); NearestMR19 = or(N22MR19,N16MR19);
                           for overlapj19=1:L19  %overlapj=overlap checking j                              
                                        if NearestMR19(round(abs(Y19(overlapj19))+0.5),round((X19(overlapj19))+0.5))==1
                                          [Ny19,Nx19]= removeoverlap(overlapj19,X19,Y19,L19,CNT19,NearestMR19);
                                          X19(overlapj19) = Nx19; Y19(overlapj19) = Ny19;
                                        end
                           end
                      [X19,Y19] = deformation(TGM19,Dmax19,X19,Y19,XX19,YY19,VVV19,k1,dEimdx,dEimdy,GM,wd19,win19,wex19,wdamp,m1,L19,bp,SVMF,SD,Nv,NearestMR19);

                      N3MR22 = roipoly(I, X3+0.5, -(Y3)+0.5);N16MR22 = roipoly(I, X16+0.5, -(Y16)+0.5);N19MR22 = roipoly(I, X19+0.5, -(Y19)+0.5);ROI22 = or(N3MR22,N16MR22);
                      NearestMR22 = or(ROI22,N19MR22);


                           for overlapj22=1:L22  %overlapj=overlap checking j 
                               %overlapj22
                                        if NearestMR22(round(abs(Y22(overlapj22))+0.5),round((X22(overlapj22))+0.5))==1
                                          [Ny22,Nx22]= removeoverlap(overlapj22,X22,Y22,L22,CNT22,NearestMR22);
                                          X22(overlapj22) = Nx22; Y22(overlapj22) = Ny22;
                                        end
                           end
                     [X22,Y22] = deformation(TGM22,Dmax22,X22,Y22,XX22,YY22,VVV22,k1,dEimdx,dEimdy,GM,wd22,win22,wex22,wdamp,m1,L22,bp,SVMF,SD,Nv,NearestMR22);
                           
%                 figure(2)
%                 imagesc(I), colormap(gray(256));
%                 hold on
%                 plot(X3+0.5,-(Y3)+0.5,'b',X16+0.5,-(Y16)+0.5,'r',X19+0.5,-(Y19)+0.5,'g',X22+0.5,-(Y22)+0.5,'y')
%                 pause
        end  % deformation loop
                %isd
end % Filter sigma loop
%end
X3=X3+0.5;Y3=Y3+0.5;X16=X16+0.5;Y16=Y16+0.5; X19=X19+0.5;Y19=Y19+0.5;X22=X22+0.5;Y22=Y22+0.5;
Y3 = -(Y3); Y16 = -(Y16);Y19 = -(Y19);Y22 = -(Y22);
% complete boundary plot
CX3=zeros(1,L3+1);CX3(1,1:L3)=X3;CX3(1,L3+1)=X3(1,1);
CY3=zeros(1,L3+1);CY3(1,1:L3)=Y3;CY3(1,L3+1)=Y3(1,1);
LL3 = length(MX3);
CMX3=zeros(1,LL3+1);CMX3(1,1:LL3)=MX3;CMX3(1,LL3+1)=MX3(1,1);
CMY3=zeros(1,LL3+1);CMY3(1,1:LL3)=MY3;CMY3(1,LL3+1)=MY3(1,1);
CXXX3=zeros(1,L3+1);CXXX3(1,1:L3)=XXX3;CXXX3(1,L3+1)=XXX3(1,1);
CYYY3=zeros(1,L3+1);CYYY3(1,1:L3)=YYY3;CYYY3(1,L3+1)=YYY3(1,1);

% figure
% imagesc(I), colormap(gray(256));
% title('With Refinement')
% hold on
% plot(CX3,CY3,'r',...
%     'LineWidth',1)
% hold on
% plot(CMX3,CMY3,'g',...
%     'LineWidth',1)
% % hold on
% % plot(X11,Y11,'b',...
% %     'LineWidth',1)
% hold on
% plot(X16,Y16,'c',...
%     'LineWidth',1)
% %axis off
% 
% figure
% imagesc(I), colormap(gray(256));
% title('Without Refinement')
% hold on
% plot(CXXX3,CYYY3,'r',...
%     'LineWidth',1)
% hold on
% plot(CMX3,CMY3,'g',...
%     'LineWidth',1)
% hold on
% plot(MX16,MY16,'g',...
%     'LineWidth',1)
% hold on
% plot(MX19,MY19,'g',...
%     'LineWidth',1)
% hold on
% plot(MX3,MY3,'g',...
%     'LineWidth',1)
% hold on
% plot(X3,Y3,'b',...
%     'LineWidth',1)
% hold on
% plot(X16,Y16,'c',...
%     'LineWidth',1)
% hold on
% plot(X17,Y17,'y',...
%     'LineWidth',1)
% hold on
% plot(X19,Y19,'m',...
%     'LineWidth',1)
% hold on
% plot(X25,Y25,'b',...
%     'LineWidth',1)
%axis off

%plot(CXXX,CYYY,'r',MX,MY,'g')
% 
% figure
% imagesc(GM), colormap(gray(256));
% hold on
% plot(CX,CY,'b',CXXX,CYYY,'r')


% DSC calculation for refinement
%yy = spline(X3,Y3,xx)
%CIX3 = interp1(1:length(CX3),CX3,1:0.2:length(CX3),'pchip');
%CIY3 = interp1(1:length(CY3),CY3,1:0.2:length(CY3),'pchip');
%A = roipoly(I,CIX3,r2);
A = poly2mask(CX3, CY3, 256, 256);  % Automask
%max(A(:))
%muscle_roi = inpolygon(X3,Y3,X3,Y3)
%CIMX3 = interp1(1:length(CMX3),CMX3,1:0.2:length(CMX3),'pchip');CIMY3 = interp1(1:length(CMY3),CMY3,1:0.2:length(CMY3),'pchip');
M = poly2mask(CMX3, CMY3, 256, 256);
AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
DSCRefinement3 = ((2*NAM)/(NA+NM))

A = poly2mask(X16, Y16, 256, 256); % Automask
M = poly2mask(MX16, MY16, 256, 256);
AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
DSCRefinement16 = ((2*NAM)/(NA+NM))

A = poly2mask(X19, Y19, 256, 256); % Automask
M = poly2mask(MX19, MY19, 256, 256);
AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
DSCRefinement19 = ((2*NAM)/(NA+NM))

A = poly2mask(X22, Y22, 256, 256); % Automask
M = poly2mask(MX22, MY22, 256, 256);
AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
DSCRefinement22 = ((2*NAM)/(NA+NM))

% DSC calculation for label transfer
A = poly2mask(XXX3, YYY3, 256, 256); % Automask
M = poly2mask(MX3, MY3, 256, 256);
AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
DSC2Label_T3 = ((2*NAM)/(NA+NM))

A = poly2mask(XXX16, YYY16, 256, 256); % Automask
M = poly2mask(MX16, MY16, 256, 256);
AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
DSC2Label_T16 = ((2*NAM)/(NA+NM))

A = poly2mask(XXX19, YYY19, 256, 256); % Automask
M = poly2mask(MX19, MY19, 256, 256);
AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
DSC2Label_T19 = ((2*NAM)/(NA+NM))

A = poly2mask(XXX22, YYY22, 256, 256); % Automask
M = poly2mask(MX22, MY22, 256, 256);
AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
DSC2Label_T22 = ((2*NAM)/(NA+NM))





    


