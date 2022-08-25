% Grouping Deformation Algorithm, saving for ROI(14:25); taking from
% structure
%3 force act; internal(local curvature), Gaussian potential, external
%dynamic distance, 1st search along local curvature then along reverse
%local curvature, then calculate acceleration , velocity, position change,
%If current poins is boundary then no change, If don't get boundary point
%within search then acts two force(internal+external). filter applied in
%four stages varing sigma. Remove overlap in each iteration for each muscle

clc
clear all 
close all
global R_muscle_data
% Loading Nobel image volume and ground truth
tic
load('Patient 1 original MRI.mat') %contains CT, ct-info
load('R-muscleDataP2toP1.mat')
[M,N,O] = size(CT);


% for mm = 1:25
% 	for ll = 1:2
% 		for ss = 1:O
% 			R_muscle_data{mm}{ll}{ss}.x = []; % R_muscle_data= refined muscle data storing cell array
% 			R_muscle_data{mm}{ll}{ss}.y = [];
% 			R_muscle_data{mm}{ll}{ss}.saved = 0;
% 		end
% 	end
% end
% save R-muscleDataP41toP25 R_muscle_data %R-muscleDataP2toP1= refined muscle data from patient2 to patient1

% Loading transfered labels
load('T-muscleDataP2toP1.mat')
for ll = 1:2
    for ss = 14:25
        II = zeros(1024,1024);
        II(:,:) = CT(:,:,ss);  % Nobel Image
        I = imresize(II,0.25); %256,256
        
        L3=0;L11=0;L16=0;L19=0;L22=0;L24=0;
        X3 = T_muscle_data{3}{ll}{ss}.x;LL3 = length(X3);Y3= T_muscle_data{3}{ll}{ss}.y;X11 = T_muscle_data{11}{ll}{ss}.x;LL11 = length(X11);Y11= T_muscle_data{11}{ll}{ss}.y;
        X16 = T_muscle_data{16}{ll}{ss}.x;LL16 = length(X16);Y16= T_muscle_data{16}{ll}{ss}.y;X19 = T_muscle_data{19}{ll}{ss}.x;LL19 = length(X19);Y19= T_muscle_data{19}{ll}{ss}.y;
        X22 = T_muscle_data{22}{ll}{ss}.x;LL22 = length(X22);Y22= T_muscle_data{22}{ll}{ss}.y;X24 = T_muscle_data{24}{ll}{ss}.x;LL24 = length(X24);Y24= T_muscle_data{24}{ll}{ss}.y;
%         figure(1)
%         imagesc(I), colormap(gray(256));
%         hold on
%         plot(X3,Y3,'b',MX3,MY3,'g',X16,Y16,'r',MX16,MY16,'g',X19,Y19,'c',MX19,MY19,'g',X22,Y22,'y',MX22,MY22,'g')

        %%%%%%%%%%%%%%%%%%
        if LL3>=1
        %RESAMPLING muscle3
        x = 1:length(X3);v = X3;xq = 1:0.5:length(X3);X3=[];
        X3 = interp1(x,v,xq,'spline');                       
        x = 1:length(Y3);v = Y3;xq = 1:0.5:length(Y3);Y3=[];
        Y3 = interp1(x,v,xq,'spline');       
        %Centroid Calculation muscle3
        ROIM3 = roipoly(I, X3, Y3); %ROIM =ROI Musk
        CENS3 = regionprops(ROIM3,'Centroid'); %CENS =Centroid Structure
        CNTy3 =  CENS3(1,1).Centroid(2)-0.5;CNTy3 = -(CNTy3);
        CNT3 =[CENS3(1,1).Centroid(1)-0.5,CNTy3]; %CNT=Centroid
        
        X3=X3-0.5;Y3=Y3-0.5;Y3 = -(Y3);L3 = length(X3);  % Transfering contour into cartesian coordinate
        end
        %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%
        if LL11>=1
        x = 1:length(X11);v = X11;xq = 1:0.5:length(X11);X11=[];
        X11 = interp1(x,v,xq,'spline');                       %RESAMPLING muscle11
        x = 1:length(Y11);v = Y11;xq = 1:0.5:length(Y11);Y11=[];
        Y11 = interp1(x,v,xq,'spline');
        
        ROIM11 = roipoly(I, X11, Y11); CENS11 = regionprops(ROIM11,'Centroid'); 
        CNTy11 =  CENS11(1,1).Centroid(2)-0.5;CNTy11 = -(CNTy11);CNT11 =[CENS11(1,1).Centroid(1)-0.5,CNTy11];
        X11=X11-0.5;Y11=Y11-0.5;Y11 = -(Y11);L11 = length(X11);        
        end
        %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%
        if LL16>=1
        x = 1:length(X16);v = X16;xq = 1:0.5:length(X16);X16=[];
        X16 = interp1(x,v,xq,'spline');                       %RESAMPLING muscle16
        x = 1:length(Y16);v = Y16;xq = 1:0.5:length(Y16);Y16=[];
        Y16 = interp1(x,v,xq,'spline');
        
        ROIM16 = roipoly(I, X16, Y16); CENS16 = regionprops(ROIM16,'Centroid'); CNTy16 =  CENS16(1,1).Centroid(2)-0.5;
        CNTy16 = -(CNTy16);CNT16 =[CENS16(1,1).Centroid(1)-0.5,CNTy16];
        X16=X16-0.5;Y16=Y16-0.5;Y16 = -(Y16);L16 = length(X16);
        end
        %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%
        if LL19>=1
        x = 1:length(X19);v = X19;xq = 1:0.5:length(X19);X19=[];
        X19 = interp1(x,v,xq,'spline');                       %RESAMPLING muscle19
        x = 1:length(Y19);v = Y19;xq = 1:0.5:length(Y19);Y19=[];
        Y19 = interp1(x,v,xq,'spline');
        
        ROIM19 = roipoly(I, X19, Y19); CENS19 = regionprops(ROIM19,'Centroid'); CNTy19 =  CENS19(1,1).Centroid(2)-0.5;
        CNTy19 = -(CNTy19);CNT19 =[CENS19(1,1).Centroid(1)-0.5,CNTy19];
        X19=X19-0.5;Y19=Y19-0.5;Y19 = -(Y19);L19 = length(X19);
        end
        %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%
        if LL22>=1
        x = 1:length(X22);v = X22;xq = 1:0.5:length(X22);X22=[];
        X22 = interp1(x,v,xq,'spline');                       %RESAMPLING muscle22
        x = 1:length(Y22);v = Y22;xq = 1:0.5:length(Y22);Y22=[];
        Y22 = interp1(x,v,xq,'spline');
        
%         imagesc(I); colormap(gray(256));
%         hold on 
%         plot(X22,Y22,'g')
        
        ROIM22 = roipoly(I, X22, Y22);
        CENS22 = regionprops(ROIM22,'Centroid'); 
        CNTy22 =  CENS22(1,1).Centroid(2)-0.5;
        CNTy22 = -(CNTy22);CNT22 =[CENS22(1,1).Centroid(1)-0.5,CNTy22];
        X22=X22-0.5;Y22=Y22-0.5;Y22 = -(Y22);L22 = length(X22);
        end
        %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%
        if LL24>=1
        x = 1:length(X24);v = X24;xq = 1:0.5:length(X24);X24=[];
        X24 = interp1(x,v,xq,'spline');                       %RESAMPLING muscle3
        x = 1:length(Y24);v = Y24;xq = 1:0.5:length(Y24);Y24=[];
        Y24 = interp1(x,v,xq,'spline');
        
        ROIM24 = roipoly(I, X24, Y24); CENS24 = regionprops(ROIM24,'Centroid'); CNTy24 =  CENS24(1,1).Centroid(2)-0.5;
        CNTy24 = -(CNTy24);CNT24 =[CENS24(1,1).Centroid(1)-0.5,CNTy24];
        X24=X24-0.5;Y24=Y24-0.5;Y24 = -(Y24);L24 = length(X24);
        end
        %%%%%%%%%%%%%%%%%%%%%%
%         % saving initial positions
%         XXX3 = X3;
%         YYY3 = Y3; 
%         XXX16 = X16;XXX19 = X19;XXX22 = X22;YYY16 = Y16;YYY19 = Y19;YYY22 = Y22;              
        % ROI= or(ROIM3,ROIM16);
        % ROI1= or(ROI,ROIM19);
        % figure
        % imagesc(ROI1), colormap(gray(256));
%         figure(2)
%         imagesc(I), colormap(gray(256));
%         hold on
%         plot(X3+0.5,-(Y3)+0.5,'b',X16+0.5,-(Y16)+0.5,'r',X19+0.5,-(Y19)+0.5,'g',X22+0.5,-(Y22)+0.5,'y')

        sigma=9;TGM3 =630;TGM11 =2;TGM16 =10;TGM19 =5;TGM22 =7;TGM24 =1; Dmax3 = 20;Dmax11 = 26;Dmax16 = 20;Dmax19 = 20;Dmax22 = 20;Dmax24 = 20;
        win3 = 100;win11 = 1000;win16 = 1000;win19 = 1000;win22 = 1000;win24 = 1000; % internal force waiting
        wex3 = 50;wex11 = 500;wex16 = 50;wex19 = 50;wex22 = 50;wex24 = 50;wd3 = 5.8;wd11 = 17.8;wd16 = 3.8;wd19 = 17.8;wd22 = 3.8;wd24 = 17.8; %im=20;%wdamp=1.6;
                %TSDmax = 0;
        for isd=1:4 %isd=iteration starnard diviation;  filter varing loop
           k1 = [0 0 0 0 0 0 -.5 1 -.5 0 0 0 0 0 0];
           win3 = win3/10000;win11 = win11/10000;win16 = win16/10000;win19 = win19/10000;win22 = win22/10000;win24 = win24/10000; wex3 = wex3/10000;wex11 = wex11/10000;wex16 = wex16/10000;wex19 = wex19/10000;wex22 = wex22/10000;wex24 = wex24/10000;
           wd3 = wd3/3;wd11 = wd11/3;wd16 = wd16/3;wd19 = wd19/3;wd22 = wd22/3;wd24 = wd24/3; wdamp = -.8;  
           Dmax3 = ceil(Dmax3/3);Dmax11 = ceil(Dmax11/3);Dmax16 = ceil(Dmax16/3);Dmax19 = ceil(Dmax19/3);Dmax22 = ceil(Dmax22/3);Dmax24 = ceil(Dmax24/3); % distance threshold
           TGM3 =3*(isd)*TGM3;TGM11 =3*(isd)*TGM11;TGM16 =3*(isd)*TGM16;TGM19 =3*(isd)*TGM19;TGM22 =3*(isd)*TGM22;TGM24 =3*(isd)*TGM24; % TGM= threshold gradient magnitude
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
          XX3 = zeros(1,L3);XX11 = zeros(1,L11);XX16 = zeros(1,L16);XX19 = zeros(1,L19);XX22 = zeros(1,L22);XX24 = zeros(1,L24); % saving deformed contour x values for all vertices
          YY3 = zeros(1,L3);YY11 = zeros(1,L11);YY16 = zeros(1,L16);YY19 = zeros(1,L19);YY22 = zeros(1,L22);YY24 = zeros(1,L24);
          VVV3 = zeros(L3,2);VVV11 = zeros(L11,2);VVV16 = zeros(L16,2);VVV19 = zeros(L19,2);VVV22 = zeros(L22,2);VVV24 = zeros(L24,2); % initial velocities       
          %SV = zeros(1,Dmax); % SV=search vector        

               for i=1:10 % deformation iteration
                            %i
                            R3 = roipoly(I, X3+0.5, -(Y3)+0.5);R11 = roipoly(I, X11+0.5, -(Y11)+0.5);R16 = roipoly(I, X16+0.5, -(Y16)+0.5);
                            R19 = roipoly(I, X19+0.5, -(Y19)+0.5);R22 = roipoly(I, X22+0.5, -(Y22)+0.5);R24 = roipoly(I, X24+0.5, -(Y24)+0.5);
                            if L3>=1
                               NearestMR3 = R22; % NearestMR3 = Neatrest muscle region for m3 creating binary mask, 1 inside                     
                               % Checking m3 overlap
                                   for overlapj3=1:L3  %overlapj=overlap checking j                              
                                            if NearestMR3(round(abs(Y3(overlapj3))+0.5),round(X3(overlapj3)+0.5))==1
                                              [Ny3,Nx3]= removeoverlap(overlapj3,X3,Y3,L3,CNT3,NearestMR3);
                                              X3(overlapj3) = Nx3; Y3(overlapj3) = Ny3;
                                            end
                                   end
                               [X3,Y3] = deformation(TGM3,Dmax3,X3,Y3,XX3,YY3,VVV3,k1,dEimdx,dEimdy,GM,wd3,win3,wex3,wdamp,m1,L3,bp,SVMF,SD,Nv,NearestMR3);
                            end
                            if L11>=1                     
                               ROI11 = or(R16,R19);NearestMR11 = or(ROI11,R24);
                                   for overlapj11=1:L11  %overlapj=overlap checking j                              
                                            if NearestMR11(round(abs(Y11(overlapj11))+0.5),round(X11(overlapj11)+0.5))==1
                                              [Ny11,Nx11]= removeoverlap(overlapj11,X11,Y11,L11,CNT11,NearestMR11);
                                              X11(overlapj11) = Nx11; Y11(overlapj11) = Ny11;
                                            end
                                   end
                               [X11,Y11] = deformation(TGM11,Dmax11,X11,Y11,XX11,YY11,VVV11,k1,dEimdx,dEimdy,GM,wd11,win11,wex11,wdamp,m1,L11,bp,SVMF,SD,Nv,NearestMR11);
                            end
                            if L16>=1
                              ROI161 = or(R22,R19); ROI162 = or(R11,R24);NearestMR16 = or(ROI161,ROI162);
%                                figure(3)
%                                imagesc(I), colormap(gray(256));
%                                hold on
%                                plot(X3+0.5,-(Y3)+0.5,'b',X16+0.5,-(Y16)+0.5,'r',X19+0.5,-(Y19)+0.5,'g',X22+0.5,-(Y22)+0.5,'y')
                                  %pause
                                  for overlapj16=1:L16  %overlapj=overlap checking j                              
                                            if NearestMR16(round(abs(Y16(overlapj16))+0.5),round((X16(overlapj16))+0.5))==1
                                              [Ny16,Nx16]= removeoverlap(overlapj16,X16,Y16,L16,CNT16,NearestMR16);
                                              X16(overlapj16) = Nx16; Y16(overlapj16) = Ny16;
                                            end
                                  end
%                                figure(4)
%                                imagesc(I), colormap(gray(256));
%                                hold on
%                                plot(X3+0.5,-(Y3)+0.5,'b',X16+0.5,-(Y16)+0.5,'r',X19+0.5,-(Y19)+0.5,'g',X22+0.5,-(Y22)+0.5,'y')
                                   %pause
                               [X16,Y16] = deformation(TGM16,Dmax16,X16,Y16,XX16,YY16,VVV16,k1,dEimdx,dEimdy,GM,wd16,win16,wex16,wdamp,m1,L16,bp,SVMF,SD,Nv,NearestMR16);
                            end
                            
                            if L19>=1
                               ROI191 = or(R22,R16); ROI192 = or(R11,R24);NearestMR19 = or(ROI191,ROI192);
                                   for overlapj19=1:L19  %overlapj=overlap checking j                              
                                                if NearestMR19(round(abs(Y19(overlapj19))+0.5),round((X19(overlapj19))+0.5))==1
                                                  [Ny19,Nx19]= removeoverlap(overlapj19,X19,Y19,L19,CNT19,NearestMR19);
                                                  X19(overlapj19) = Nx19; Y19(overlapj19) = Ny19;
                                                end
                                   end
                              [X19,Y19] = deformation(TGM19,Dmax19,X19,Y19,XX19,YY19,VVV19,k1,dEimdx,dEimdy,GM,wd19,win19,wex19,wdamp,m1,L19,bp,SVMF,SD,Nv,NearestMR19);
                            end
                            if L22>=1
                              ROI221 = or(R19,R16); ROI222 = or(R3,R24);NearestMR22 = or(ROI221,ROI222);
                                   for overlapj22=1:L22  %overlapj=overlap checking j 
                                       %overlapj22
                                                if NearestMR22(round(abs(Y22(overlapj22))+0.5),round((X22(overlapj22))+0.5))==1
                                                  [Ny22,Nx22]= removeoverlap(overlapj22,X22,Y22,L22,CNT22,NearestMR22);
                                                  X22(overlapj22) = Nx22; Y22(overlapj22) = Ny22;
                                                end
                                   end
                             [X22,Y22] = deformation(TGM22,Dmax22,X22,Y22,XX22,YY22,VVV22,k1,dEimdx,dEimdy,GM,wd22,win22,wex22,wdamp,m1,L22,bp,SVMF,SD,Nv,NearestMR22);
                            end
                             if L24>=1
                                ROI241 = or(R22,R16); ROI242 = or(R11,R19);NearestMR24 = or(ROI241,ROI242);
                                   for overlapj24=1:L24  %overlapj=overlap checking j                              
                                            if NearestMR24(round(abs(Y24(overlapj24))+0.5),round(X24(overlapj24)+0.5))==1
                                              [Ny24,Nx24]= removeoverlap(overlapj24,X24,Y24,L24,CNT24,NearestMR24);
                                              X24(overlapj24) = Nx24; Y24(overlapj24) = Ny24;
                                            end
                                   end
                               [X24,Y24] = deformation(TGM24,Dmax24,X24,Y24,XX24,YY24,VVV24,k1,dEimdx,dEimdy,GM,wd24,win24,wex24,wdamp,m1,L24,bp,SVMF,SD,Nv,NearestMR24);
                            end

        %                 figure(2)
        %                 imagesc(I), colormap(gray(256));
        %                 hold on
        %                 plot(X3+0.5,-(Y3)+0.5,'b',X16+0.5,-(Y16)+0.5,'r',X19+0.5,-(Y19)+0.5,'g',X22+0.5,-(Y22)+0.5,'y')
        %                 pause
                end  % deformation loop
                        %isd
        end % Filter sigma loop
        % Saving Refined data
          if L3>=1
              X3=X3+0.5;Y3=Y3+0.5;Y3 = -(Y3);R_muscle_data{3}{ll}{ss}.x = X3; R_muscle_data{3}{ll}{ss}.y = Y3;R_muscle_data{3}{ll}{ss}.saved = 1;             
          end
          if L11>=1
              X11=X11+0.5;Y11=Y11+0.5;Y11 = -(Y11);R_muscle_data{11}{ll}{ss}.x = X11; R_muscle_data{11}{ll}{ss}.y = Y11;R_muscle_data{11}{ll}{ss}.saved = 1;             
          end
          if L16>=1
              X16=X16+0.5;Y16=Y16+0.5; Y16 = -(Y16);R_muscle_data{16}{ll}{ss}.x = X16; R_muscle_data{16}{ll}{ss}.y = Y16;R_muscle_data{16}{ll}{ss}.saved = 1;
          end
          if L19>=1
              X19=X19+0.5;Y19=Y19+0.5; Y19 = -(Y19);R_muscle_data{19}{ll}{ss}.x = X19; R_muscle_data{19}{ll}{ss}.y = Y19;R_muscle_data{19}{ll}{ss}.saved = 1;
          end
          if L22>=1
              X22=X22+0.5;Y22=Y22+0.5;Y22 = -(Y22);R_muscle_data{22}{ll}{ss}.x = X22; R_muscle_data{22}{ll}{ss}.y = Y22;R_muscle_data{22}{ll}{ss}.saved = 1;
          end
        if L24>=1
             X24=X24+0.5;Y24=Y24+0.5;Y24 = -(Y24); R_muscle_data{24}{ll}{ss}.x = X24; R_muscle_data{24}{ll}{ss}.y = Y24;R_muscle_data{24}{ll}{ss}.saved = 1;
        end
              %L3=0;L11=0;L16=0;L19=0;L22=0;L24=0;
              %ss
    end
    
end
save R-muscleDataP2toP1 R_muscle_data
toc
% % complete boundary plot
% CX3=zeros(1,L3+1);CX3(1,1:L3)=X3;CX3(1,L3+1)=X3(1,1);
% CY3=zeros(1,L3+1);CY3(1,1:L3)=Y3;CY3(1,L3+1)=Y3(1,1);
% LL3 = length(MX3);
% CMX3=zeros(1,LL3+1);CMX3(1,1:LL3)=MX3;CMX3(1,LL3+1)=MX3(1,1);
% CMY3=zeros(1,LL3+1);CMY3(1,1:LL3)=MY3;CMY3(1,LL3+1)=MY3(1,1);
% CXXX3=zeros(1,L3+1);CXXX3(1,1:L3)=XXX3;CXXX3(1,L3+1)=XXX3(1,1);
% CYYY3=zeros(1,L3+1);CYYY3(1,1:L3)=YYY3;CYYY3(1,L3+1)=YYY3(1,1);
% 


% % figure
% % imagesc(I), colormap(gray(256));
% % title('With Refinement')
% % hold on
% % plot(CX3,CY3,'r',...
% %     'LineWidth',1)
% % hold on
% % plot(CMX3,CMY3,'g',...
% %     'LineWidth',1)
% % % hold on
% % % plot(X11,Y11,'b',...
% % %     'LineWidth',1)
% % hold on
% % plot(X16,Y16,'c',...
% %     'LineWidth',1)
% % %axis off
% % 
% % figure
% % imagesc(I), colormap(gray(256));
% % title('Without Refinement')
% % hold on
% % plot(CXXX3,CYYY3,'r',...
% %     'LineWidth',1)
% % hold on
% % plot(CMX3,CMY3,'g',...
% %     'LineWidth',1)
% % hold on
% % plot(MX16,MY16,'g',...
% %     'LineWidth',1)
% % hold on
% % plot(MX19,MY19,'g',...
% %     'LineWidth',1)
% % hold on
% % plot(MX3,MY3,'g',...
% %     'LineWidth',1)
% % hold on
% % plot(X3,Y3,'b',...
% %     'LineWidth',1)
% % hold on
% % plot(X16,Y16,'c',...
% %     'LineWidth',1)
% % hold on
% % plot(X17,Y17,'y',...
% %     'LineWidth',1)
% % hold on
% % plot(X19,Y19,'m',...
% %     'LineWidth',1)
% % hold on
% % plot(X25,Y25,'b',...
% %     'LineWidth',1)
% %axis off
% 
% %plot(CXXX,CYYY,'r',MX,MY,'g')
% % 
% % figure
% % imagesc(GM), colormap(gray(256));
% % hold on
% % plot(CX,CY,'b',CXXX,CYYY,'r')
% 
% 
% % DSC calculation for refinement
% %yy = spline(X3,Y3,xx)
% %CIX3 = interp1(1:length(CX3),CX3,1:0.2:length(CX3),'pchip');
% %CIY3 = interp1(1:length(CY3),CY3,1:0.2:length(CY3),'pchip');
% %A = roipoly(I,CIX3,r2);
% A = poly2mask(CX3, CY3, 256, 256);  % Automask
% %max(A(:))
% %muscle_roi = inpolygon(X3,Y3,X3,Y3)
% %CIMX3 = interp1(1:length(CMX3),CMX3,1:0.2:length(CMX3),'pchip');CIMY3 = interp1(1:length(CMY3),CMY3,1:0.2:length(CMY3),'pchip');
% M = poly2mask(CMX3, CMY3, 256, 256);
% AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
% DSCRefinement3 = ((2*NAM)/(NA+NM))
% 
% A = poly2mask(X16, Y16, 256, 256); % Automask
% M = poly2mask(MX16, MY16, 256, 256);
% AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
% DSCRefinement16 = ((2*NAM)/(NA+NM))
% 
% A = poly2mask(X19, Y19, 256, 256); % Automask
% M = poly2mask(MX19, MY19, 256, 256);
% AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
% DSCRefinement19 = ((2*NAM)/(NA+NM))
% 
% A = poly2mask(X22, Y22, 256, 256); % Automask
% M = poly2mask(MX22, MY22, 256, 256);
% AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
% DSCRefinement22 = ((2*NAM)/(NA+NM))
% 
% % DSC calculation for label transfer
% A = poly2mask(XXX3, YYY3, 256, 256); % Automask
% M = poly2mask(MX3, MY3, 256, 256);
% AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
% DSC2Label_T3 = ((2*NAM)/(NA+NM))
% 
% A = poly2mask(XXX16, YYY16, 256, 256); % Automask
% M = poly2mask(MX16, MY16, 256, 256);
% AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
% DSC2Label_T16 = ((2*NAM)/(NA+NM))
% 
% A = poly2mask(XXX19, YYY19, 256, 256); % Automask
% M = poly2mask(MX19, MY19, 256, 256);
% AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
% DSC2Label_T19 = ((2*NAM)/(NA+NM))
% 
% A = poly2mask(XXX22, YYY22, 256, 256); % Automask
% M = poly2mask(MX22, MY22, 256, 256);
% AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
% DSC2Label_T22 = ((2*NAM)/(NA+NM))





    


