% Deformation Algorithm
%3 force act; internal(local curvature), Gaussian potential, external
%dynamic distance, 1st search along local curvature then along reverse
%local curvature, then ccalculate acceleration , velocity, position change,
%If current poins is boundary then no change, If don't get boundary point
%within search then acts two force(internal+external). filter applied in
%four stages varing sigma




clc
clear all 
close all
% load('TMD-m11sRf11.mat')
% X11 = xnnnnn; % X11=11 no muscle
% Y11= ynnnnn; 
% load('TMD-m16sLf10.mat')
% X16 = xnnnnn; % X11=11 no muscle
% Y16= ynnnnn; 

muscle = [3 11 16 22];
frames = [16:20; 10:14; 9:13; 8:12];
muscleSide = {'L','R'};

for musS = 1:2
    for mus = 4
    mus;
            for frm = 1:5
                frm;
                load(['TMD-m',num2str(muscle(mus)),'s',muscleSide{musS},'f',num2str(frames(mus,frm)),'.mat'])  % contains xnnnn, ynnnn,xnnnnn,ynnnnn
                load('P1.mat')
                load('p2 muscle data.mat') %p2=P1, contains manual contour , can be used for comparison
                MX = muscle_data{MN}{SNP}{SNIRI}.x; % manual data for comparison; 1=left
                MY = muscle_data{MN}{SNP}{SNIRI}.y;
        %         MX16 = muscle_data{16}{SNP}{SNIRI}.x; % manual data for comparison; 1=left
        %         MY16 = muscle_data{16}{SNP}{SNIRI}.y;
                I = zeros(256,256);
                I(:,:) = V(:,:,SNIRI);
                %SNIRI
                % figure
                % imagesc(I), colormap(gray(256));
                % [D1,D2] = size(I);
                % [xx,yy] = meshgrid([1:D2]-120,[1:D1]-123);
                %a = min(I)

                % Contour x y values retriving
                % XXXX = muscle_data{5}{1}{2}.x;  % 5th muscle 1= left side, 1=slice x data
                % YYYY = muscle_data{5}{1}{2}.y;  % Y= y coordinates values for all vertices
                % LL = length(XXXX);
                % X = XXXX(1:LL-1);
                % Y = YYYY(1:LL-1);
                  X = xnnnnn; % without local registration "xnnnn"
                  Y= ynnnnn;  % without local registration "xnnnnn"
                % imagesc(I), colormap(gray(256));
                % hold on
                % plot(X,Y,'b',xxx,yyy,'r')

                %%%%%%%%%%%%%%%%%%
                x = 1:length(X);
                v = X;
                xq = 1:0.5:length(X);
                X=[];
                X = interp1(x,v,xq,'spline');                       %RESAMPLING

                x = 1:length(Y);
                v = Y;
                xq = 1:0.5:length(Y);
                Y=[];
                Y = interp1(x,v,xq,'spline');
                %%%%%%%%%%%%%%%%%%%%%%

                XXX = X;
                YYY = Y; % saving initial positions
                %Centroid Calculation
                ROIM = roipoly(I, X, Y); %ROIM =ROI Musk
                CENS = regionprops(ROIM,'Centroid'); %CENS =Centroid Structure
                CNTy =  CENS.Centroid(2)-0.5;
                CNTy = -(CNTy);
                CNT =[CENS.Centroid(1)-0.5,CNTy]; %CNT=Centroid

                X=X-0.5;
                Y=Y-0.5;
                Y = -(Y);
                L = length(X);
                % X(19)
                % Y(19)
                % Centroid Calculation
                
                
                
                %image features for distance force
%                 G1=fspecial('gauss',[40 40],0.1);
%                 CI1(:,:)=filter2(G1,I);[dx1,dy1] = gradient(CI1); dx1 = dx1.^2; dy1 = dy1.^2; GM = dx1+dy1;

                % XX = zeros(L,2); % XX= initial contour matrix
                % for i = 1:L
                %     XX(i,1) = X(i);
                %     XX(i,2) = Y(i);
                % end
                sigma=9;TGM =.1;TGM2 = 18000; Dmax = 20;
                win = 100000; % internal force waiting
                wex = 50;wd = 17.8; %im=20;%wdamp=1.6;
                %TSDmax = 0;
                for isd=1:4 %isd=iteration starnard diviation;  filter varing loop
                    isd;
                        k1 = [0 0 0 0 0 0 -.5 1 -.5 0 0 0 0 0 0];
                        win = win/10000; wex = wex/10000;wd = wd/3; wdamp = -.8;  
                        Dmax = ceil(Dmax/3); % distance threshold
                        TGM =3*(isd)*TGM;%TGM2 =3*(isd)*TGM2; % TGM= threshold gradient magnitude
                        sigma=sigma/3;
                        %im=im/2;

                        % Calculating potential energy function/distribution
                        G=fspecial('gauss',[40 40],sigma);
                        CI(:,:)=filter2(G,I);[dCIdx,dCIdy] = gradient(CI); dx = dCIdx.^2; dy = dCIdy.^2; GM = dx+dy;
                        GM1 = dCIdx+dCIdy;[dCIdx1,dCIdy1] = gradient(GM1);GM2 = dCIdx1+dCIdy1;
                        %GM = CI;
                        % figure
                        % imagesc(GM), colormap(gray(256));
                        Eim = -(GM2); % - means valley
                        [dEimdx,dEimdy] = gradient(Eim);
                        dEimdx = -1*dEimdx; % external force field
                        dEimdy = -1*dEimdy; 
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Starting Deformation process
                        XX = zeros(1,L); % deformed contour x values for all vertices
                        YY = zeros(1,L);
                        VVV = zeros(L,2); % initial velocities       
                        %SV = zeros(1,Dmax); % SV=search vector        
                        bp = [0 0]; % bp=boundary point
                        Nv = [0;0]; % Nv=normal vector
                        m1 = 1; % m1= mass for all vertices
                        SVMF=0; % stop vertex moving flag
                        SD = 0;
                        for i=1:10 % deformation iteration
                             %i
%                              ROIM = roipoly(I, X+0.5, -Y+0.5); %ROIM =ROI Musk
%                              CENS = regionprops(ROIM,'Centroid'); %CENS =Centroid Structure
%                              CNT =[CENS.Centroid(1),CENS.Centroid(2)]; %CNT=Centroid
                                    for j=1:L % vertex selection loop
                                           %j
                                            P1 = [X(j),Y(j)];  % position vectors
                                                if j==L
                                                    P2 = [X(1),Y(1)];   % Initial nth, 1st, 2nd vertices
                                                else
                                                    P2 = [X(j+1),Y(j+1)];    
                                                end

                                                if j==1
                                                    Pn = [X(L),Y(L)];
                                                else
                                                    Pn = [X(j-1),Y(j-1)]; 
                                                end
                                            v1 = [VVV(j,1);VVV(j,2)]; 

                                            % calculating Internal forces
                                            d1 = P2-P1;dn=P1-Pn;
                                            % Resampling 
                                    %         ldes = 4; lmin = 0.5*ldes; lmax = 1.5*ldes;
                                    %         
                                    %         if norm(d1)<lmin
                                    %             
                                    %         elseif norm(d1)>lmin 
                                    %             
                                    %         end
                                    %         
                                            ud1 = d1/(norm(d1)); udn = dn/(norm(dn));
                                            ut1 = (ud1+udn)/(norm(ud1+udn)); %unit tangent
                                            ur1 = [0 1;-1 0]*ut1'; % column vector ; unit normal
                                            c1 = ud1-udn; un1 = c1/(norm(c1)); % un1=unit normal; local curvature
                                            run1 = -(un1);% rotated unit normal
                                           ucd1 = (CNT-P1)/(norm(CNT-P1)); rucd1 = -(ucd1);%rucd1= riverse unit centroid direction
                                            %k1 = [0 0 0 0 0 -.0005 -.0005 .01 -.0005 -.0005 0 0 0 0 0];
                                            %dp = [c1(1)*ur1(1) c1(2)*ur1(2)];
                                            Finl1 = conv((dot(c1,ur1')),k1);
                                            %Finl1 = conv(dp,k1);
                                            %Fin1 = (norm(Finl1))*ur1; % column vector
                                            Fin1 = Finl1(8)*ur1; % column vector
                                            %Fin1 = (dot(c1,ur1))*ur1;
                                            % calculating external force on a vertex
                                    %         fx = interp2(xx,yy,dEimdx,X(1),Y(1));
                                    %         fy = interp2(xx,yy,dEimdy,X(1),Y(1));
                                            fx = interp2(dEimdx,P1(1),abs(P1(2)));
                                            fy = interp2(dEimdy,P1(1),abs(P1(2)));
                                            Fim1 = [fx,fy];
                                            %  z1 = fx+(fy*i);
                                            %  r1 = abs(z1);
                                            %  theta1 = angle(z1);
                                            Fex1 = (dot(Fim1',ur1))*ur1; % column vector

                                            %calculating Dynamic distance force
                                            theta1 = atan2d(un1(2),un1(1)); % inward searching
                                            theta2 = atan2d(run1(2),run1(1)); % outward searching
%                                             theta1 = atan2d(ucd1(2),ucd1(1)); % inward searching
%                                             theta2 = atan2d(rucd1(2),rucd1(1)); % outward searching
                                            ii = ceil(P1(2)); % row
                                            jj = ceil(P1(1));
                                            cp = [jj ii]; %cp=contour point

                                            CGM = GM(abs(ii),jj); % CGM=current GM
                                                 if CGM>TGM %&& CGM<TGM2  % checking current point as boundary point
                                                    bp = [jj ii]; % bp=boundary point
                                                    Nv = [0;0];SD=norm(cp-bp); % Nv= normar vector
                                                    SVMF=1;  %stop vertex moving flag 
                                                 end

                                             if norm(bp)==0 % Starting search

                                                            if theta1>=-5 && theta1<=5 % along curvature or normal search start
                                                                    for iii= 1:Dmax
                                                                         jj=jj+1;
                                                                         if jj>256 || jj<1 % for object near to image boundary check
                                                                             SD = 0; Nv = [0;0];break;
                                                                         end
                                                                        CGM = GM(abs(ii),jj);
                                                                        if CGM>TGM %&& CGM<TGM2
                                                                            bp = [jj ii]; % bp=boundary point
                                                                            Nv = un1'; SD=norm(cp-bp); % Nv= normar vector                                    
                                                                            break;
                                                                        end
                                                                    end
                                                            elseif theta1>=6 && theta1<=84
                                                                   for iii= 1:Dmax
                                                                         jj=jj+1;ii=ii+1;
                                                                             if jj>256 || jj<1
                                                                             SD = 0; Nv = [0;0];break;
                                                                             elseif abs(ii)>256 || abs(ii)<1
                                                                                SD = 0; Nv = [0;0]; break;
                                                                              end
                                                                         CGM = GM(abs(ii),jj);
                                                                        if CGM>TGM %&& CGM<TGM2
                                                                            bp = [jj ii]; Nv = un1';SD=norm(cp-bp); 
                                                                            break;
                                                                        end
                                                                   end
                                                            elseif theta1>=85 && theta1<=95
                                                                   for iii= 1:Dmax
                                                                         ii=ii+1;
                                                                             if abs(ii)>256 || abs(ii)<1
                                                                                SD = 0; Nv = [0;0]; break;
                                                                              end
                                                                         CGM = GM(abs(ii),jj);
                                                                        if CGM>TGM %&& CGM<TGM2
                                                                            bp = [jj ii]; Nv = un1';SD=norm(cp-bp); 
                                                                            break;
                                                                        end
                                                                   end
                                                             elseif theta1>=96 && theta1<=174
                                                                   for iii= 1:Dmax
                                                                         jj=jj-1;ii=ii+1;
                                                                           if jj>256 || jj<1
                                                                             SD = 0; Nv = [0;0];break;
                                                                          elseif abs(ii)>256 || abs(ii)<1
                                                                                 SD = 0; Nv = [0;0];break;
                                                                           end
                                                                         CGM = GM(abs(ii),jj);
                                                                        if CGM>TGM %&& CGM<TGM2
                                                                            bp = [jj ii]; Nv = un1';SD=norm(cp-bp); 
                                                                            break;
                                                                        end
                                                                   end  
                                                             elseif theta1>=175 && theta1<=-175
                                                                   for iii= 1:Dmax
                                                                         jj=jj-1;
                                                                           if jj>256 || jj<1
                                                                            SD = 0; Nv = [0;0]; break;                                                  
                                                                           end
                                                                         CGM = GM(abs(ii),jj);
                                                                        if CGM>TGM %&& CGM<TGM2
                                                                            bp = [jj ii]; Nv = un1'; SD=norm(cp-bp);
                                                                            break;
                                                                        end
                                                                   end
                                                             elseif theta1>=-174 && theta1<=-96
                                                                   for iii= 1:Dmax
                                                                         jj=jj-1; ii=ii-1;
                                                                           if jj>256 || jj<1
                                                                             SD = 0; Nv = [0;0];break;
                                                                           elseif abs(ii)>256 || abs(ii)<1
                                                                                SD = 0; Nv = [0;0]; break;
                                                                           end
                                                                         CGM = GM(abs(ii),jj);
                                                                        if CGM>TGM %&& CGM<TGM2
                                                                            bp = [jj ii];Nv = un1';SD=norm(cp-bp);
                                                                            break;
                                                                        end
                                                                   end
                                                             elseif theta1>=-95 && theta1<=-85
                                                                   for iii= 1:Dmax
                                                                        ii=ii-1;                                                  
                                                                             if abs(ii)>256 || abs(ii)<1
                                                                                SD = 0; Nv = [0;0]; break;
                                                                             end
                                                                        CGM = GM(abs(ii),jj);
                                                                        if CGM>TGM %&& CGM<TGM2
                                                                            bp = [jj ii];Nv = un1';SD=norm(cp-bp);
                                                                            break;
                                                                        end
                                                                   end
                                                             elseif theta1>=-84 && theta1<=-6
                                                                   for iii= 1:Dmax
                                                                        ii=ii-1;jj=jj+1;
                                                                          if jj>256 || jj<1
                                                                             SD = 0; Nv = [0;0];break;
                                                                          elseif abs(ii)>256 || abs(ii)<1
                                                                                SD = 0; Nv = [0;0]; break;
                                                                          end
                                                                        CGM = GM(abs(ii),jj);
                                                                        if CGM>TGM %&& CGM<TGM2
                                                                            bp = [jj ii];Nv = un1';SD=norm(cp-bp);
                                                                            break;
                                                                        end
                                                                   end
                                                                   %%%%%% along reverse curvature or normal Search start
                                                              elseif theta2>=-5 && theta2<=5
                                                                    for iii= 1:Dmax
                                                                         jj=jj+1;
                                                                           if jj>256 || jj<1
                                                                            SD = 0; Nv = [0;0]; break;                                                    
                                                                           end
                                                                        CGM = GM(abs(ii),jj); 
                                                                        if CGM>TGM %&& CGM<TGM2
                                                                            bp = [jj ii]; % bp=boundary point
                                                                            Nv = run1'; SD=norm(cp-bp); % Nv= normar vector
                                                                            break;
                                                                        end
                                                                    end
                                                            elseif theta2>=6 && theta2<=84
                                                                   for iii= 1:Dmax
                                                                         jj=jj+1;ii=ii+1;
                                                                           if jj>256 || jj<1
                                                                             SD = 0; Nv = [0;0];break;
                                                                           elseif abs(ii)>256 || abs(ii)<1
                                                                            SD = 0; Nv = [0;0]; break;
                                                                           end
                                                                         CGM = GM(abs(ii),jj);
                                                                        if CGM>TGM %&& CGM<TGM2
                                                                            bp = [jj ii]; Nv = run1';SD=norm(cp-bp); 
                                                                            break;
                                                                        end
                                                                   end
                                                            elseif theta2>=85 && theta2<=95
                                                                   for iii= 1:Dmax
                                                                         ii=ii+1;                                                  
                                                                           if abs(ii)>256 || abs(ii)<1
                                                                             SD = 0; Nv = [0;0];break;
                                                                           end
                                                                         CGM = GM(abs(ii),jj);
                                                                        if CGM>TGM %&& CGM<TGM2
                                                                            bp = [jj ii]; Nv = run1';SD=norm(cp-bp); 
                                                                            break;
                                                                        end
                                                                   end
                                                             elseif theta2>=96 && theta2<=174
                                                                   for iii= 1:Dmax
                                                                         jj=jj-1;ii=ii+1;
                                                                           if jj>256 || jj<1
                                                                             SD = 0; Nv = [0;0];break;
                                                                           elseif abs(ii)>256 || abs(ii)<1
                                                                             SD = 0; Nv = [0;0];break;
                                                                           end
                                                                         CGM = GM(abs(ii),jj);
                                                                        if CGM>TGM %&& CGM<TGM2
                                                                            bp = [jj ii]; Nv = run1';SD=norm(cp-bp); 
                                                                            break;
                                                                        end
                                                                   end  
                                                             elseif theta2>=175 && theta2<=-175
                                                                   for iii= 1:Dmax
                                                                         jj=jj-1;
                                                                           if jj>256 || jj<1
                                                                             SD = 0; Nv = [0;0];break;                                                     
                                                                           end
                                                                         CGM = GM(abs(ii),jj); 
                                                                        if CGM>TGM %&& CGM<TGM2
                                                                            bp = [jj ii]; Nv = run1';SD=norm(cp-bp); 
                                                                            break;
                                                                        end
                                                                   end
                                                             elseif theta2>=-174 && theta2<=-96
                                                                   for iii= 1:Dmax
                                                                         jj=jj-1; ii=ii-1;
                                                                           if jj>256 || jj<1
                                                                             SD = 0; Nv = [0;0];break;
                                                                           elseif abs(ii)>256 || abs(ii)<1
                                                                             SD = 0; Nv = [0;0];break;
                                                                           end
                                                                         CGM = GM(abs(ii),jj);
                                                                        if CGM>TGM %&& CGM<TGM2
                                                                            bp = [jj ii];Nv = run1';SD=norm(cp-bp);
                                                                            break;
                                                                        end
                                                                   end
                                                             elseif theta2>=-95 && theta2<=-85
                                                                   for iii= 1:Dmax
                                                                        ii=ii-1;
                                                                             if abs(ii)>256 || abs(ii)<1
                                                                                 SD = 0; Nv = [0;0];
                                                                                 break;
                                                                              end
                                                                        CGM = GM(abs(ii),jj);
                                                                        if CGM>TGM %&& CGM<TGM2
                                                                            bp = [jj ii];Nv = run1';SD=norm(cp-bp);
                                                                            break;
                                                                        end
                                                                   end
                                                             elseif theta2>=-84 && theta2<=-6
                                                                   for iii= 1:Dmax
                                                                        ii=ii-1;jj=jj+1;
                                                                          if jj>256 || jj<1
                                                                              SD = 0; Nv = [0;0];
                                                                             break;
                                                                          elseif abs(ii)>256 || abs(ii)<1
                                                                              SD = 0; Nv = [0;0];
                                                                              break;
                                                                          end
                                                                        CGM = GM(abs(ii),jj);
                                                                        if CGM>TGM %&& CGM<TGM2
                                                                            bp = [jj ii];Nv = run1'; SD=norm(cp-bp);
                                                                            break;
                                                                        end
                                                                   end
                                                             else
                                                                SD = 0; Nv = [0;0];
                                                            end 

                                             end  
                                              %%%%%%%%%%%%%%%%

                                             if SVMF==1
                                                XX(1,j)=P1(1,1);YY(1,j)=P1(1,2);
                                                VVV(j,1) = 0;VVV(j,2) = 0;bp = [0 0]; SVMF=0;
                                             else                      
                                                Fd = wd*(SD/Dmax)*Nv; % dynamic distance force                      
                                                %F1 = (Fin1*win)+(wex*Fex1)+Fd;
                                                F1 = (Fin1*win)+(wex*Fex1)+(wdamp*v1)+Fd;
                                                %F1 = (wex*Fex1)+(wdamp*v1);
                                                a1 = (1/m1)*F1; v1 = v1+a1;
                                                VVV(j,1) = v1(1,1); % needs to save v1
                                                VVV(j,2) = v1(2,1);
                                                P1 = P1+v1';                      
                                                XX(1,j)=P1(1,1);YY(1,j)=P1(1,2);bp = [0 0];%SD = 0; Nv = [0;0];
                                             end            
                                    end
                        X = XX; Y = YY;
                %         figure(1)
                %         imagesc(GM), colormap(gray(256));
                %         hold on
                %         plot(X+0.5,-(Y)+0.5,'b',XXX,YYY,'r')
                %         pause
                        end
                end
                X=X+0.5;Y=Y+0.5;
                Y = -(Y);
                

                % complete boundary plot
                CX=zeros(1,L+1);CX(1,1:L)=X;CX(1,L+1)=X(1,1);
                CY=zeros(1,L+1);CY(1,1:L)=Y;CY(1,L+1)=Y(1,1);
                LL = length(MX);
                CMX=zeros(1,LL+1);CMX(1,1:LL)=MX;CMX(1,LL+1)=MX(1,1);
                CMY=zeros(1,LL+1);CMY(1,1:LL)=MY;CMY(1,LL+1)=MY(1,1);
                CXXX=zeros(1,L+1);CXXX(1,1:L)=XXX;CXXX(1,L+1)=XXX(1,1);
                CYYY=zeros(1,L+1);CYYY(1,1:L)=YYY;CYYY(1,L+1)=YYY(1,1);

%                 figure
%                 imagesc(I), colormap(gray(256));
%                 title('With Refinement')
%                 hold on
%                 plot(CX,CY,'r',...
%                     'LineWidth',1)
%                 hold on
%                 plot(CMX,CMY,'g',...
%                     'LineWidth',1)
%                 % hold on
%                 % plot(X11,Y11,'b',...
%                 %     'LineWidth',1)
%                 % hold on
%                 % plot(X16,Y16,'c',...
%                 %     'LineWidth',1)
%                 axis off
% 
%                 figure
%                 imagesc(I), colormap(gray(256));
%                 title('Without Refinement')
%                 hold on
%                 plot(CXXX,CYYY,'r',...
%                     'LineWidth',1)
%                 hold on
%                 plot(CMX,CMY,'g',...
%                     'LineWidth',1)
%                 % hold on
%                 % plot(MX16,MY16,'g',...
%                 %     'LineWidth',1)
%                 % hold on
%                 % plot(X11,Y11,'b',...
%                 %     'LineWidth',1)
%                 % hold on
%                 % plot(X16,Y16,'c',...
%                 %     'LineWidth',1)
%                 axis off

                %plot(CXXX,CYYY,'r',MX,MY,'g')
                % 
                % figure
                % imagesc(GM), colormap(gray(256));
                % hold on
                % plot(CX,CY,'b',CXXX,CYYY,'r')


                % DSC calculation for refinement
                A = poly2mask(X, Y, 256, 256); % Automask
                M = poly2mask(MX, MY, 256, 256);
                AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
                DSCRefinement(mus,(musS-1)*5+frm) = ((2*NAM)/(NA+NM));

                % DSC calculation for label transfer
                A = poly2mask(CXXX, CYYY, 256, 256); % Automask
                M = poly2mask(MX, MY, 256, 256);
                AiM = and(A,M);A = double(A);M = double(M);AiM = double(AiM);NA= sum(A(:)); NM= sum(M(:));NAM=sum(AiM(:));
                DSC2Label_T(mus,(musS-1)*5+frm) = ((2*NAM)/(NA+NM));
            end
    end
end
DSCRefinement
 DSCRefinement-DSC2Label_T
 MRefinement = mean(DSCRefinement,2)

 
MDSC2Label_T = mean(mean(DSC2Label_T))

