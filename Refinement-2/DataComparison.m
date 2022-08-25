clc
clear all 
close all
load('Patient 1 original MRI.mat')
load('Patient 1 muscle data.mat')
load('R-muscleDataP2toP1.mat')
load('T-muscleDataP2toP1.mat')
[m,n,z] = size(CT);
II = zeros(1024,1024);

for i = 14:25
    
        II(:,:) = CT(:,:,i);
        I = imresize(II,0.25);
        x3LG = muscle_data{3}{1}{i}.x;y3LG = muscle_data{3}{1}{i}.y;x3RG = muscle_data{3}{2}{i}.x;y3RG = muscle_data{3}{2}{i}.y;
        x3LT = T_muscle_data{3}{1}{i}.x;y3LT = T_muscle_data{3}{1}{i}.y;x3RT = T_muscle_data{3}{2}{i}.x;y3RT = T_muscle_data{3}{2}{i}.y;
        x3LR = R_muscle_data{3}{1}{i}.x;y3LR = R_muscle_data{3}{1}{i}.y;x3RR = R_muscle_data{3}{2}{i}.x;y3RR = R_muscle_data{3}{2}{i}.y;
        
        x11LG = muscle_data{11}{1}{i}.x;y11LG = muscle_data{11}{1}{i}.y;x11RG = muscle_data{11}{2}{i}.x;y11RG = muscle_data{11}{2}{i}.y;
        x11LT = T_muscle_data{11}{1}{i}.x;y11LT = T_muscle_data{11}{1}{i}.y;x11RT = T_muscle_data{11}{2}{i}.x;y11RT = T_muscle_data{11}{2}{i}.y;
        x11LR = R_muscle_data{11}{1}{i}.x;y11LR = R_muscle_data{11}{1}{i}.y;x11RR = R_muscle_data{11}{2}{i}.x;y11RR = R_muscle_data{11}{2}{i}.y;
        
        x16LG = muscle_data{16}{1}{i}.x;y16LG = muscle_data{16}{1}{i}.y;x16RG = muscle_data{16}{2}{i}.x;y16RG = muscle_data{16}{2}{i}.y;
        x16LT = T_muscle_data{16}{1}{i}.x;y16LT = T_muscle_data{16}{1}{i}.y;x16RT = T_muscle_data{16}{2}{i}.x;y16RT = T_muscle_data{16}{2}{i}.y;
        x16LR = R_muscle_data{16}{1}{i}.x;y16LR = R_muscle_data{16}{1}{i}.y;x16RR = R_muscle_data{16}{2}{i}.x;y16RR = R_muscle_data{16}{2}{i}.y;
        
        x19LG = muscle_data{19}{1}{i}.x;y19LG = muscle_data{19}{1}{i}.y;x19RG = muscle_data{19}{2}{i}.x;y19RG = muscle_data{19}{2}{i}.y;
        x19LT = T_muscle_data{19}{1}{i}.x;y19LT = T_muscle_data{19}{1}{i}.y;x19RT = T_muscle_data{19}{2}{i}.x;y19RT = T_muscle_data{19}{2}{i}.y;
        x19LR = R_muscle_data{19}{1}{i}.x;y19LR = R_muscle_data{19}{1}{i}.y;x19RR = R_muscle_data{19}{2}{i}.x;y19RR = R_muscle_data{19}{2}{i}.y;
        
        x22LG = muscle_data{22}{1}{i}.x;y22LG = muscle_data{22}{1}{i}.y;x22RG = muscle_data{22}{2}{i}.x;y22RG = muscle_data{22}{2}{i}.y;
        x22LT = T_muscle_data{22}{1}{i}.x;y22LT = T_muscle_data{22}{1}{i}.y;x22RT = T_muscle_data{22}{2}{i}.x;y22RT = T_muscle_data{22}{2}{i}.y;
        x22LR = R_muscle_data{22}{1}{i}.x;y22LR = R_muscle_data{22}{1}{i}.y;x22RR = R_muscle_data{22}{2}{i}.x;y22RR = R_muscle_data{22}{2}{i}.y;
        
        x24LG = muscle_data{24}{1}{i}.x;y24LG = muscle_data{24}{1}{i}.y;x24RG = muscle_data{24}{2}{i}.x;y24RG = muscle_data{24}{2}{i}.y;
        x24LT = T_muscle_data{24}{1}{i}.x;y24LT = T_muscle_data{24}{1}{i}.y;x24RT = T_muscle_data{24}{2}{i}.x;y24RT = T_muscle_data{24}{2}{i}.y;
        x24LR = R_muscle_data{24}{1}{i}.x;y24LR = R_muscle_data{24}{1}{i}.y;x24RR = R_muscle_data{24}{2}{i}.x;y24RR = R_muscle_data{24}{2}{i}.y;
        
       % % complete boundary plot
        LLG3 = length(x3LG);LLT3 = length(x3LT);LLR3 = length(x3LR); LRG3 = length(x3RG);LRT3 = length(x3RT);LRR3 = length(x3RR);
        LLG11 = length(x11LG);LLT11 = length(x11LT);LLR11 = length(x11LR); LRG11 = length(x11RG);LRT11 = length(x11RT);LRR11 = length(x11RR);
        LLG16 = length(x16LG);LLT16 = length(x16LT);LLR16 = length(x16LR); LRG16 = length(x16RG);LRT16 = length(x16RT);LRR16 = length(x16RR);
        LLG22 = length(x22LG);LLT22 = length(x22LT);LLR22 = length(x22LR); LRG22 = length(x22RG);LRT22 = length(x22RT);LRR22 = length(x22RR);
          
        Cx3LG=zeros(1,LLG3+1);Cx3LG(1,1:LLG3)=x3LG;Cx3LG(1,LLG3+1)=x3LG(1,1); Cx3LT=zeros(1,LLT3+1);Cx3LT(1,1:LLT3)=x3LT;Cx3LT(1,LLT3+1)=x3LT(1,1); Cx3LR=zeros(1,LLR3+1);Cx3LR(1,1:LLR3)=x3LR;Cx3LR(1,LLR3+1)=x3LR(1,1);
        Cy3LG=zeros(1,LLG3+1);Cy3LG(1,1:LLG3)=y3LG;Cy3LG(1,LLG3+1)=y3LG(1,1); Cy3LT=zeros(1,LLT3+1);Cy3LT(1,1:LLT3)=y3LT;Cy3LT(1,LLT3+1)=y3LT(1,1); Cy3LR=zeros(1,LLR3+1);Cy3LR(1,1:LLR3)=y3LR;Cy3LR(1,LLR3+1)=y3LR(1,1);

        Cx3RG=zeros(1,LRG3+1);Cx3RG(1,1:LRG3)=x3RG;Cx3RG(1,LRG3+1)=x3RG(1,1); Cx3RT=zeros(1,LRT3+1);Cx3RT(1,1:LRT3)=x3RT;Cx3RT(1,LRT3+1)=x3RT(1,1); Cx3RR=zeros(1,LRR3+1);Cx3RR(1,1:LRR3)=x3RR;Cx3RR(1,LRR3+1)=x3RR(1,1);
        Cy3RG=zeros(1,LRG3+1);Cy3RG(1,1:LRG3)=y3RG;Cy3RG(1,LRG3+1)=y3RG(1,1); Cy3RT=zeros(1,LRT3+1);Cy3RT(1,1:LRT3)=y3RT;Cy3RT(1,LRT3+1)=y3RT(1,1); Cy3RR=zeros(1,LRR3+1);Cy3RR(1,1:LRR3)=y3RR;Cy3RR(1,LRR3+1)=y3RR(1,1);
        
        Cx11LG=zeros(1,LLG11+1);Cx11LG(1,1:LLG11)=x11LG;Cx11LG(1,LLG11+1)=x11LG(1,1); Cx11LT=zeros(1,LLT11+1);Cx11LT(1,1:LLT11)=x11LT;Cx11LT(1,LLT11+1)=x11LT(1,1); Cx11LR=zeros(1,LLR11+1);Cx11LR(1,1:LLR11)=x11LR;Cx11LR(1,LLR11+1)=x11LR(1,1);
        Cy11LG=zeros(1,LLG11+1);Cy11LG(1,1:LLG11)=y11LG;Cy11LG(1,LLG11+1)=y11LG(1,1); Cy11LT=zeros(1,LLT11+1);Cy11LT(1,1:LLT11)=y11LT;Cy11LT(1,LLT11+1)=y11LT(1,1); Cy11LR=zeros(1,LLR11+1);Cy11LR(1,1:LLR11)=y11LR;Cy11LR(1,LLR11+1)=y11LR(1,1);

        Cx11RG=zeros(1,LRG11+1);Cx11RG(1,1:LRG11)=x11RG;Cx11RG(1,LRG11+1)=x11RG(1,1); Cx11RT=zeros(1,LRT11+1);Cx11RT(1,1:LRT11)=x11RT;Cx11RT(1,LRT11+1)=x11RT(1,1); Cx11RR=zeros(1,LRR11+1);Cx11RR(1,1:LRR11)=x11RR;Cx11RR(1,LRR11+1)=x11RR(1,1);
        Cy11RG=zeros(1,LRG11+1);Cy11RG(1,1:LRG11)=y11RG;Cy11RG(1,LRG11+1)=y11RG(1,1); Cy11RT=zeros(1,LRT11+1);Cy11RT(1,1:LRT11)=y11RT;Cy11RT(1,LRT11+1)=y11RT(1,1); Cy11RR=zeros(1,LRR11+1);Cy11RR(1,1:LRR11)=y11RR;Cy11RR(1,LRR11+1)=y11RR(1,1);
        
        Cx16LG=zeros(1,LLG16+1);Cx16LG(1,1:LLG16)=x16LG;Cx16LG(1,LLG16+1)=x16LG(1,1); Cx16LT=zeros(1,LLT16+1);Cx16LT(1,1:LLT16)=x16LT;Cx16LT(1,LLT16+1)=x16LT(1,1); Cx16LR=zeros(1,LLR16+1);Cx16LR(1,1:LLR16)=x16LR;Cx16LR(1,LLR16+1)=x16LR(1,1);
        Cy16LG=zeros(1,LLG16+1);Cy16LG(1,1:LLG16)=y16LG;Cy16LG(1,LLG16+1)=y16LG(1,1); Cy16LT=zeros(1,LLT16+1);Cy16LT(1,1:LLT16)=y16LT;Cy16LT(1,LLT16+1)=y16LT(1,1); Cy16LR=zeros(1,LLR16+1);Cy16LR(1,1:LLR16)=y16LR;Cy16LR(1,LLR16+1)=y16LR(1,1);

        Cx16RG=zeros(1,LRG16+1);Cx16RG(1,1:LRG16)=x16RG;Cx16RG(1,LRG16+1)=x16RG(1,1); Cx16RT=zeros(1,LRT16+1);Cx16RT(1,1:LRT16)=x16RT;Cx16RT(1,LRT16+1)=x16RT(1,1); Cx16RR=zeros(1,LRR16+1);Cx16RR(1,1:LRR16)=x16RR;Cx16RR(1,LRR16+1)=x16RR(1,1);
        Cy16RG=zeros(1,LRG16+1);Cy16RG(1,1:LRG16)=y16RG;Cy16RG(1,LRG16+1)=y16RG(1,1); Cy16RT=zeros(1,LRT16+1);Cy16RT(1,1:LRT16)=y16RT;Cy16RT(1,LRT16+1)=y16RT(1,1); Cy16RR=zeros(1,LRR16+1);Cy16RR(1,1:LRR16)=y16RR;Cy16RR(1,LRR16+1)=y16RR(1,1);
        
        Cx22LG=zeros(1,LLG22+1);Cx22LG(1,1:LLG22)=x22LG;Cx22LG(1,LLG22+1)=x22LG(1,1); Cx22LT=zeros(1,LLT22+1);Cx22LT(1,1:LLT22)=x22LT;Cx22LT(1,LLT22+1)=x22LT(1,1); Cx22LR=zeros(1,LLR22+1);Cx22LR(1,1:LLR22)=x22LR;Cx22LR(1,LLR22+1)=x22LR(1,1);
        Cy22LG=zeros(1,LLG22+1);Cy22LG(1,1:LLG22)=y22LG;Cy22LG(1,LLG22+1)=y22LG(1,1); Cy22LT=zeros(1,LLT22+1);Cy22LT(1,1:LLT22)=y22LT;Cy22LT(1,LLT22+1)=y22LT(1,1); Cy22LR=zeros(1,LLR22+1);Cy22LR(1,1:LLR22)=y22LR;Cy22LR(1,LLR22+1)=y22LR(1,1);

        Cx22RG=zeros(1,LRG22+1);Cx22RG(1,1:LRG22)=x22RG;Cx22RG(1,LRG22+1)=x22RG(1,1); Cx22RT=zeros(1,LRT22+1);Cx22RT(1,1:LRT22)=x22RT;Cx22RT(1,LRT22+1)=x22RT(1,1); Cx22RR=zeros(1,LRR22+1);Cx22RR(1,1:LRR22)=x22RR;Cx22RR(1,LRR22+1)=x22RR(1,1);
        Cy22RG=zeros(1,LRG22+1);Cy22RG(1,1:LRG22)=y22RG;Cy22RG(1,LRG22+1)=y22RG(1,1); Cy22RT=zeros(1,LRT22+1);Cy22RT(1,1:LRT22)=y22RT;Cy22RT(1,LRT22+1)=y22RT(1,1); Cy22RR=zeros(1,LRR22+1);Cy22RR(1,1:LRR22)=y22RR;Cy22RR(1,LRR22+1)=y22RR(1,1);
        
        
        % Plotting
        imagesc(I), colormap(gray(256));
        hold on
        
        % For refinement and result check up
%         plot(Cx3LT,Cy3LT,'g',Cx3LR,Cy3LR,'r',Cx3RT,Cy3RT,'g',Cx3RR,Cy3RR,'r','LineWidth',1.5)
%         %plot(Cx11LT,Cy11LT,'g',Cx11LR,Cy11LR,'b',Cx11RT,Cy11RT,'g',Cx11RR,Cy11RR,'b','LineWidth',1.5)
%         plot(Cx16LT,Cy16LT,'g',Cx16LR,Cy16LR,'m',Cx16RT,Cy16RT,'g',Cx16RR,Cy16RR,'m','LineWidth',1.5)
%         %plot(x19LG,y19LG,'r',x19LT,y19LT,'g',x19LR,y19LR,'b',x19RG,y19RG,'r',x19RT,y19RT,'g',x19RR,y19RR,'b')
%         plot(Cx22LT,Cy22LT,'g',Cx22LR,Cy22LR,'y',Cx22RT,Cy22RT,'g',Cx22RR,Cy22RR,'y','LineWidth',1.5)
%        % plot(x24LG,y24LG,'r',x24LT,y24LT,'g',x24LR,y24LR,'b',x24RG,y24RG,'r',x24RT,y24RT,'g',x24RR,y24RR,'b')
%        axis off
%         pause
        
        
        % For Auto segmentation
        plot(Cx3LG,Cy3LG,'g',Cx3LR,Cy3LR,'r',Cx3RG,Cy3RG,'g',Cx3RR,Cy3RR,'r','LineWidth',1.5)
        plot(Cx11LG,Cy11LG,'g',Cx11LR,Cy11LR,'b',Cx11RG,Cy11RG,'g',Cx11RR,Cy11RR,'b','LineWidth',1.5)
        plot(Cx16LG,Cy16LG,'g',Cx16LR,Cy16LR,'m',Cx16RG,Cy16RG,'g',Cx16RR,Cy16RR,'m','LineWidth',1.5)
        %plot(x19LG,y19LG,'g',x19LT,y19LT,'g',x19LR,y19LR,'b',x19RG,y19RG,'r',x19RT,y19RT,'g',x19RR,y19RR,'b')
        plot(Cx22LG,Cy22LG,'g',Cx22LR,Cy22LR,'y',Cx22RG,Cy22RG,'g',Cx22RR,Cy22RR,'y','LineWidth',1.5)
       % plot(x24LG,y24LG,'g',x24LT,y24LT,'g',x24LR,y24LR,'b',x24RG,y24RG,'r',x24RT,y24RT,'g',x24RR,y24RR,'b')
       axis off
        pause
end

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
% for i=1:z
%   % Retyriving x y data for left side muscles  
% xl1 = muscle_data{1}{1}{i}.x;
% yl1 = muscle_data{1}{1}{i}.y;
% xl2 = muscle_data{2}{1}{i}.x;
% yl2 = muscle_data{2}{1}{i}.y;
% % xl3 = muscle_data{3}{1}{z}.x;
% % yl3 = muscle_data{3}{1}{z}.y;
% % xl4 = muscle_data{4}{1}{z}.x;
% % yl4 = muscle_data{4}{1}{z}.y;
% xl5 = muscle_data{5}{1}{i}.x; % centre
% yl5 = muscle_data{5}{1}{i}.y;
% % xl7 = muscle_data{7}{1}{z}.x;
% % yl7 = muscle_data{7}{1}{z}.y;
% xrl1 = xl1+256;xrl5 = xl5+256; xrl2 = xl2+256;%xrl3 = xl3+256; xrl4 = xl4+256; xrl7 = xl7+256;
% 
% 
% % Retyriving x y data for right side muscles
% 
% xr1 = muscle_data{1}{2}{i}.x;
% yr1 = muscle_data{1}{2}{i}.y;
% xr2 = muscle_data{2}{2}{i}.x;
% yr2 = muscle_data{2}{2}{i}.y;
% % xr3 = muscle_data{3}{2}{z}.x;
% % yr3 = muscle_data{3}{2}{z}.y;
% % xr4 = muscle_data{4}{2}{z}.x;
% % yr4 = muscle_data{4}{2}{z}.y;
% % xr5 = muscle_data{5}{2}{z}.x;
% % yr5 = muscle_data{5}{2}{z}.y;
% xrr1 = xr1+256; xrr2 = xr2+256;%xrr3 = xr3+256; xrr4 = xr4+256;xrr5 = xr5+256;
%     %%%%%%%%%
%     ct1=x1(:,:,i);
%     ct2=x2(:,:,i);
% %     imagesc([ct1 ct2]), colormap(gray(256));
% imagesc([ct2 ct1]), colormap(gray(256));
% hold on
% plot(xl1,yl1,'g',xl5,yl5,'r',xr1,yr1,'b',xl2,yl2,'c',xr2,yr2,'m')
% plot(xrl1,yl1,'g',xrl5,yl5,'r',xrr1,yr1,'b',xrl2,yl2,'c',xrr2,yr2,'m')
% 
%  
% % imagesc(ct2); %colormap(gray(25);
% %     
%     pause
%     
% end
% 
% 
