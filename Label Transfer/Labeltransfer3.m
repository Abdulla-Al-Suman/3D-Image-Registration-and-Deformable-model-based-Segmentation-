% 45 slices, all label transfer
clc
clear all 
close all
global MN SNIRI T_muscle_data
load('Patient 41 muscle data.mat')% contains muscle_data cell array
load('T-muscleDataP41toP25.mat') % loading transfered cell array to store remaining muscles data
%load('Patient2.mat')
%load('ISBI-ini40i.mat')
load('elastic8VP41VVP25.mat')
load('elastic18VP41VVP25.mat')
load('elastic32VP41VVP25.mat')% contains V VV
load('elasticLVP41VVP25.mat')
x1=V;[M,N,O] = size(x1);



% for mm = 1:25
% 	for ll = 1:2
% 		for ss = 1:O
% 			T_muscle_data{mm}{ll}{ss}.x = [];
% 			T_muscle_data{mm}{ll}{ss}.y = [];  % Need to active and deactive			
% 			T_muscle_data{mm}{ll}{ss}.saved = 0;
% 		end
% 	end
% end
% save T-muscleDataP41toP25 T_muscle_data 



%MN = 24; SNA = 1; FN = 12; SNP=1; % MN=muscle number, SNA=1(right)= side number atlas, FN= frame number , SNP=2(right)side number patient
muscle = [3 11 16 19 22 24]; % muscles number for segmentation

for MN = 5:6
    for SNA = 1:2
        for FN = 14:28  % 1:O for full data
            if SNA == 1
            SNP=1; %SNP=2(right)side number patient
            elseif SNA == 2
                SNP=2;
            end
            xl1 = muscle_data{muscle(MN)}{SNA}{FN}.x;  % muscle selection
            yl1 = muscle_data{muscle(MN)}{SNA}{FN}.y;
            L=length(xl1);
            if L>=1         
                    xl1 = xl1-129.5;yl1 = yl1-129.5;
                    FNIA = (FN-1)-(O/2); % FNIA = frame number in affine registration %[x,y,z] = meshgrid([0:D2-1]-D2/2,[0:D1-1]-D1/2,[0:D3-1]-D3/2);
                    zl1 = FNIA*(ones(1,L));

                    % figure
                    % imagesc(V(:,:,23)), colormap(gray(256));
                    % hold on
                    % plot(xl1+128.5,yl1+128.5,'g')

                    % data transform for affine transformation
                    xn =  0.9995*xl1+0.0414*yl1+0.2497*zl1-6.0862;
                    yn = -0.0670*xl1+0.9478*yl1+0.1092*zl1-5.7872;
                    zn = -0.0033*xl1-0.0205*yl1+1.1560*zl1-1.4480;
                    %SNI = mode(zn);% have to use 44 slice, SNI=slice nimber index
                    SNI = mean(zn);% have to use 44 slice, SNI=slice nimber index
                    SNIRI = round((O/2)+SNI+1); % SNIRI=slice number in registered image 
                    xn=xn+128;yn=yn+128; % 

                    % figure
                    % imagesc([VV(:,:,SNIRI) V(:,:,SNIRI)]), colormap(gray(256));
                    % hold on
                    % plot(xn+.5,yn+.5,'g',xn+.5+256,yn+.5,'r')

                    if SNIRI>=1 && SNIRI<=O % O= total no of slices
                            % data transform for elastic transformation p=8
                            m = MFEF8(SNIRI,:);%[-0.6354    4.5398   -1.3903  -17.0288    1.1858    0.5914   -4.6825    9.1254];
                            xnn = zeros(1,L);ynn = zeros(1,L);
                            s=2;dxdmX=0;kk= 1;
                            for u = 0:s-1
                                for v = 0:s-1
                                    syms  x y
                                    dxdmX = dxdmX+m(kk)*cos(((2*x+1)*pi*u)/(2*N)).*cos(((2*y+1)*pi*v)/(2*M));
                                    kk=kk+1;
                                end
                            end

                            dxdmY=0;
                            for u = 0:s-1
                                for v = 0:s-1
                                    syms  x y
                                    dxdmY = dxdmY+m(kk)*cos(((2*x+1)*pi*u)/(2*N)).*cos(((2*y+1)*pi*v)/(2*M));
                                    kk=kk+1;
                            %         pause
                                end
                            end

                            for i=1:L
                            syms  x y
                            % [x,y] = solve(x+m(1)*cos(((2*x+1)*pi*0)/(2*N))*cos(((2*y+1)*pi*0)/(2*M))+m(2)*cos(((2*x+1)*pi*0)/(2*N))*cos(((2*y+1)*pi*1)/(2*M))+m(3)*cos(((2*x+1)*pi*1)/(2*N))*cos(((2*y+1)*pi*0)/(2*M))+m(4)*cos(((2*x+1)*pi*1)/(2*N))*cos(((2*y+1)*pi*1)/(2*M))==xn(i),y+m(5)*cos(((2*x+1)*pi*0)/(2*N))*cos(((2*y+1)*pi*0)/(2*M))+m(6)*cos(((2*x+1)*pi*0)/(2*N))*cos(((2*y+1)*pi*1)/(2*M))+m(7)*cos(((2*x+1)*pi*1)/(2*N))*cos(((2*y+1)*pi*0)/(2*M))+m(8)*cos(((2*x+1)*pi*1)/(2*N))*cos(((2*y+1)*pi*1)/(2*M))==yn(i));
                            %[x,y] = solve(x+m(1)*cos(((2*x+1)*pi*0)/2*N)*cos(((2*y+1)*pi*0)/2*M)+m(2)*cos(((2*x+1)*pi*0)/2*N)*cos(((2*y+1)*pi*1)/2*M)==xn(i),y+m(5)*cos(((2*x+1)*pi*0)/2*N)*cos(((2*y+1)*pi*0)/2*M)+m(6)*cos(((2*x+1)*pi*0)/2*N)*cos(((2*y+1)*pi*1)/2*M)==yn(i));
                            [x,y] = solve(x+dxdmX==xn(i),y+dxdmY==yn(i));
                            xnn(i) = x;ynn(i) = y;
                            end

                            % data transform for elastic transformation p=18
                            m=MFEF18(SNIRI,:);%[1.1167    0.6659    1.6352    1.3772    0.2859    1.3901   -1.0030   -3.9360   -1.6428    1.4436    1.6569    2.3429    5.5437    6.8240    7.1798   -1.3613 -3.7167   -2.3846];
                            xnnn = zeros(1,L);ynnn = zeros(1,L);
                            s=3;dxdmX=0;kk= 1;
                            for u = 0:s-1
                                for v = 0:s-1
                                    syms  x y
                                    dxdmX = dxdmX+m(kk)*cos(((2*x+1)*pi*u)/(2*N)).*cos(((2*y+1)*pi*v)/(2*M));
                                    kk=kk+1;
                                end
                            end

                            dxdmY=0;
                            for u = 0:s-1
                                for v = 0:s-1
                                    syms  x y
                                    dxdmY = dxdmY+m(kk)*cos(((2*x+1)*pi*u)/(2*N)).*cos(((2*y+1)*pi*v)/(2*M));
                                    kk=kk+1;
                            %         pause
                                end
                            end

                            for i=1:L
                            syms  x y
                            % [x,y] = solve(x+m(1)*cos(((2*x+1)*pi*0)/(2*N))*cos(((2*y+1)*pi*0)/(2*M))+m(2)*cos(((2*x+1)*pi*0)/(2*N))*cos(((2*y+1)*pi*1)/(2*M))+m(3)*cos(((2*x+1)*pi*1)/(2*N))*cos(((2*y+1)*pi*0)/(2*M))+m(4)*cos(((2*x+1)*pi*1)/(2*N))*cos(((2*y+1)*pi*1)/(2*M))==xn(i),y+m(5)*cos(((2*x+1)*pi*0)/(2*N))*cos(((2*y+1)*pi*0)/(2*M))+m(6)*cos(((2*x+1)*pi*0)/(2*N))*cos(((2*y+1)*pi*1)/(2*M))+m(7)*cos(((2*x+1)*pi*1)/(2*N))*cos(((2*y+1)*pi*0)/(2*M))+m(8)*cos(((2*x+1)*pi*1)/(2*N))*cos(((2*y+1)*pi*1)/(2*M))==yn(i));
                            %[x,y] = solve(x+m(1)*cos(((2*x+1)*pi*0)/2*N)*cos(((2*y+1)*pi*0)/2*M)+m(2)*cos(((2*x+1)*pi*0)/2*N)*cos(((2*y+1)*pi*1)/2*M)==xn(i),y+m(5)*cos(((2*x+1)*pi*0)/2*N)*cos(((2*y+1)*pi*0)/2*M)+m(6)*cos(((2*x+1)*pi*0)/2*N)*cos(((2*y+1)*pi*1)/2*M)==yn(i));
                            [x,y] = solve(x+dxdmX==xnn(i),y+dxdmY==ynn(i));
                            xnnn(i) = x;ynnn(i) = y;
                            end

                            % data transform for elastic transformation p=32
                            m=MFEF32(SNIRI,:);%[2.8157    6.6548    3.3056    2.0782   -7.4710   -7.9816  -10.7929   -2.6952    2.2324    9.0834    2.4141    4.0951   -5.7683   -6.4903  -10.7645   -2.7367  3.5738    6.4365    4.5532    0.9127   -8.5700  -16.1042   -9.0753   -6.5879    3.6837    7.5753    4.8209    1.1522   -5.2047   -7.6519   -4.7079   -2.6983];
                            xnnnn = zeros(1,L);ynnnn = zeros(1,L);
                            s=4;dxdmX=0;kk= 1;
                            for u = 0:s-1
                                for v = 0:s-1
                                    syms  x y
                                    dxdmX = dxdmX+m(kk)*cos(((2*x+1)*pi*u)/(2*N)).*cos(((2*y+1)*pi*v)/(2*M));
                                    kk=kk+1;
                                end
                            end

                            dxdmY=0;
                            for u = 0:s-1
                                for v = 0:s-1
                                    syms  x y
                                    dxdmY = dxdmY+m(kk)*cos(((2*x+1)*pi*u)/(2*N)).*cos(((2*y+1)*pi*v)/(2*M));
                                    kk=kk+1;
                            %         pause
                                end
                            end

                            for i=1:L
                            syms  x y
                            % [x,y] = solve(x+m(1)*cos(((2*x+1)*pi*0)/(2*N))*cos(((2*y+1)*pi*0)/(2*M))+m(2)*cos(((2*x+1)*pi*0)/(2*N))*cos(((2*y+1)*pi*1)/(2*M))+m(3)*cos(((2*x+1)*pi*1)/(2*N))*cos(((2*y+1)*pi*0)/(2*M))+m(4)*cos(((2*x+1)*pi*1)/(2*N))*cos(((2*y+1)*pi*1)/(2*M))==xn(i),y+m(5)*cos(((2*x+1)*pi*0)/(2*N))*cos(((2*y+1)*pi*0)/(2*M))+m(6)*cos(((2*x+1)*pi*0)/(2*N))*cos(((2*y+1)*pi*1)/(2*M))+m(7)*cos(((2*x+1)*pi*1)/(2*N))*cos(((2*y+1)*pi*0)/(2*M))+m(8)*cos(((2*x+1)*pi*1)/(2*N))*cos(((2*y+1)*pi*1)/(2*M))==yn(i));
                            %[x,y] = solve(x+m(1)*cos(((2*x+1)*pi*0)/2*N)*cos(((2*y+1)*pi*0)/2*M)+m(2)*cos(((2*x+1)*pi*0)/2*N)*cos(((2*y+1)*pi*1)/2*M)==xn(i),y+m(5)*cos(((2*x+1)*pi*0)/2*N)*cos(((2*y+1)*pi*0)/2*M)+m(6)*cos(((2*x+1)*pi*0)/2*N)*cos(((2*y+1)*pi*1)/2*M)==yn(i));
                            [x,y] = solve(x+dxdmX==xnnn(i),y+dxdmY==ynnn(i));
                            xnnnn(i) = x;ynnnn(i) = y;
                            end

                            % data transform for Local elastic transformation 
                            xnnnnn = zeros(1,L);ynnnnn = zeros(1,L);
                            for j=1:L
                                if (xnnnn(j)>=0 && xnnnn(j)<127) && (ynnnn(j)>=0 && ynnnn(j)<127)
                                        m=MFEFL(SNIRI,:);
                                        s=3;dxdmX=0;kk= 1;
                                        for u = 0:s-1
                                            for v = 0:s-1
                                                syms  x y
                                                dxdmX = dxdmX+m(kk)*cos(((2*x+1)*pi*u)/(2*128)).*cos(((2*y+1)*pi*v)/(2*128));
                                                kk=kk+1;
                                            end
                                        end
                                        dxdmY=0;
                                        for u = 0:s-1
                                            for v = 0:s-1
                                                syms  x y
                                                dxdmY = dxdmY+m(kk)*cos(((2*x+1)*pi*u)/(2*128)).*cos(((2*y+1)*pi*v)/(2*128));
                                                kk=kk+1;
                                            end
                                        end
                                        syms  x y
                                        [x,y] = solve(x+dxdmX==xnnnn(j),y+dxdmY==ynnnn(j));
                                        xnnnnn(j) = x;ynnnnn(j) = y;

                                elseif (xnnnn(j)>=127) && (ynnnn(j)>=0 && ynnnn(j)<127)
                                        m=MFEFL(SNIRI+45,:);
                                        s=3;dxdmX=0;kk= 1;
                                        for u = 0:s-1
                                            for v = 0:s-1
                                                syms  x y
                                                dxdmX = dxdmX+m(kk)*cos(((2*x+1)*pi*u)/(2*128)).*cos(((2*y+1)*pi*v)/(2*128));
                                                kk=kk+1;
                                            end
                                        end
                                        dxdmY=0;
                                        for u = 0:s-1
                                            for v = 0:s-1
                                                syms  x y
                                                dxdmY = dxdmY+m(kk)*cos(((2*x+1)*pi*u)/(2*128)).*cos(((2*y+1)*pi*v)/(2*128));
                                                kk=kk+1;
                                            end
                                        end
                                        syms  x y
                                        [x,y] = solve(x+dxdmX==xnnnn(j)-127,y+dxdmY==ynnnn(j));
                                        xnnnnn(j) = x+127;ynnnnn(j) = y;
                                elseif (xnnnn(j)>=0 && xnnnn(j)<127) && (ynnnn(j)>=127)
                                        m=MFEFL(SNIRI+90,:);
                                        s=3;dxdmX=0;kk= 1;
                                        for u = 0:s-1
                                            for v = 0:s-1
                                                syms  x y
                                                dxdmX = dxdmX+m(kk)*cos(((2*x+1)*pi*u)/(2*128)).*cos(((2*y+1)*pi*v)/(2*128));
                                                kk=kk+1;
                                            end
                                        end
                                        dxdmY=0;
                                        for u = 0:s-1
                                            for v = 0:s-1
                                                syms  x y
                                                dxdmY = dxdmY+m(kk)*cos(((2*x+1)*pi*u)/(2*128)).*cos(((2*y+1)*pi*v)/(2*128));
                                                kk=kk+1;
                                            end
                                        end
                                        syms  x y
                                        [x,y] = solve(x+dxdmX==xnnnn(j),y+dxdmY==ynnnn(j)-127);
                                        xnnnnn(j) = x;ynnnnn(j) = y+127;

                                elseif (xnnnn(j)>=127) && (ynnnn(j)>=127)
                                        m=MFEFL(SNIRI+135,:);
                                        s=3;dxdmX=0;kk= 1;
                                        for u = 0:s-1
                                            for v = 0:s-1
                                                syms  x y
                                                dxdmX = dxdmX+m(kk)*cos(((2*x+1)*pi*u)/(2*128)).*cos(((2*y+1)*pi*v)/(2*128));
                                                kk=kk+1;
                                            end
                                        end
                                        dxdmY=0;
                                        for u = 0:s-1
                                            for v = 0:s-1
                                                syms  x y
                                                dxdmY = dxdmY+m(kk)*cos(((2*x+1)*pi*u)/(2*128)).*cos(((2*y+1)*pi*v)/(2*128));
                                                kk=kk+1;
                                            end
                                        end
                                        syms  x y
                                        [x,y] = solve(x+dxdmX==xnnnn(j)-127,y+dxdmY==ynnnn(j)-127);
                                        xnnnnn(j) = x+127;ynnnnn(j) = y+127;
                                end

                            end
                            xnnnnn=xnnnnn+1.5;ynnnnn=ynnnnn+1.5;
                            % Saving transformend data
                            T_muscle_data{muscle(MN)}{SNP}{SNIRI}.x = xnnnnn; 
                            T_muscle_data{muscle(MN)}{SNP}{SNIRI}.y = ynnnnn;
                            T_muscle_data{muscle(MN)}{SNP}{SNIRI}.saved = 1;
                            
                            save T-muscleDataP41toP25 T_muscle_data
                    end % if end; 1:45 check
            end % Muscle data checking end
            FN
        end % no of slices end
        SNA
    end
    MN
end

    % T-muscleDataP2toP1 = transform muscle data from patient2 to patient1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xnnnn=xnnnn+1.5;ynnnn=ynnnn+1.5;

%save TMD-m25sLf13 xnnnnn ynnnnn xnnnn ynnnn SNIRI MN SNP

% %%%%% Plotting transformed data with ground truth
% load('p2 muscle data.mat') %p2=P1, contains manual contour , can be used for comparison
% MX = muscle_data{MN}{SNP}{SNIRI}.x; 
% MY = muscle_data{MN}{SNP}{SNIRI}.y;
% load('ISBIelastic32.mat')% contains V VV
% figure
% imagesc([VV(:,:,SNIRI) V(:,:,SNIRI)]), colormap(gray(256));
% hold on
% plot(xnnnn,ynnnn,'g',MX,MY,'b',xnnnn+256,ynnnn,'r')
% 
% % DSC calculation
% A = poly2mask(xnnnn, ynnnn, 256, 256); % Automask
% M = poly2mask(MX, MY, 256, 256);
% AiM = and(A,M);
% A = double(A);
% M = double(M);
% AiM = double(AiM);
% NA= sum(A(:)); 
% NM= sum(M(:));
% NAM=sum(AiM(:));
% DSC = ((2*NAM)/(NA+NM))





