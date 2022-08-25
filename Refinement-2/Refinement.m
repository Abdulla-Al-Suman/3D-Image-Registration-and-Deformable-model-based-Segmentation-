clc
clear all 
close all

load('p2 muscle data.mat')   % muscle_data=cell array
load('FR881832-80808060m10304050.mat')
I = zeros(256,256);
I(:,:) = V(:,:,3);
% figure
% imagesc(I), colormap(gray(256));
% [D1,D2] = size(I);
% [xx,yy] = meshgrid([1:D2]-120,[1:D1]-123);
%a = min(I)

% Contour x y values retriving
XXXX = muscle_data{5}{1}{3}.x;  % 5th muscle 1= left side, 1=slice x data
YYYY = muscle_data{5}{1}{3}.y;  % Y= y coordinates values for all vertices
%%%%%%%%%%%%%%%%%%
% x = 1:length(X);
% v = X;
% xq = 1:0.1:length(X);
% X=[];
% X = interp1(x,v,xq,'spline');                       RESAMPLING
% 
% x = 1:length(Y);
% v = Y;
% xq = 1:0.1:length(Y);
% Y=[];
% Y = interp1(x,v,xq,'spline');
%%%%%%%%%%%%%%%%%%%%%%
LL = length(XXXX);
X = XXXX(1:LL-1);
Y = YYYY(1:LL-1);

XXX = X;
YYY = Y; % saving initial positions
% imagesc(I), colormap(gray(256));
% hold on
% plot(X,Y,'b')
X=X-0.5;
Y=Y-0.5;
Y = -(Y);
L = length(X);

% XX = zeros(L,2); % XX= initial contour matrix
% for i = 1:L
%     XX(i,1) = X(i);
%     XX(i,2) = Y(i);
% end

% Calculating potential energy function/distribution
G=fspecial('gauss',[20 20],3);
CI(:,:)=filter2(G,I);
[dCIdx,dCIdy] = gradient(CI);
dx = dCIdx.^2;
dy = dCIdy.^2;
GM = dx+dy;
% figure
% imagesc(GM), colormap(gray(256));
Eim = -(GM); % - means valley
[dEimdx,dEimdy] = gradient(Eim);
dEimdx = -1*dEimdx; % external force field
dEimdy = -1*dEimdy; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting Deformation process
XX = zeros(1,L); % deformed contour x values for all vertices
YY = zeros(1,L);
VVV = zeros(L,2); % initial velocities
Dmax = 22; % distance threshold
%SV = zeros(1,Dmax); % SV=search vector
TGM =8.5; % TGM= threshold gradient magnitude
bp = [0 0]; % bp=boundary point
Nv = [0;0]; % Nv=normal vector
m1 = 1; % m1= mass for all vertices
SVMF=0; % stop vertex moving flag
SD = 0;
for i=1:40 % deformation iteration
    
    for j=1:L % vertex selection loop

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

            d1 = P2-P1;
            dn=P1-Pn;
            % Resampling 
    %         ldes = 4; lmin = 0.5*ldes; lmax = 1.5*ldes;
    %         
    %         if norm(d1)<lmin
    %             
    %         elseif norm(d1)>lmin 
    %             
    %         end
    %         
            ud1 = d1/(norm(d1));
            udn = dn/(norm(dn));
            ut1 = (ud1+udn)/(norm(ud1+udn));
            ur1 = [0 1;-1 0]*ut1'; % column vector
            c1 = ud1-udn;
            un1 = c1/(norm(c1)); % un1=unit normal
            run1 = -(un1);% rotated unit normal

            %k1 = [0 0 0 0 0 -.0005 -.0005 .01 -.0005 -.0005 0 0 0 0 0];
            k1 = [0 0 0 0 0 0 -.5 1 -.5 0 0 0 0 0 0];
            %dp = [c1(1)*ur1(1) c1(2)*ur1(2)];

            Finl1 = conv((dot(c1,ur1')),k1);
            %Finl1 = conv(dp,k1);
            %Fin1 = (norm(Finl1))*ur1; % column vector
            Finl1(8);
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
            ii = ceil(P1(2)); % row
            jj = ceil(P1(1));
            cp = [jj ii]; %cp=contour point

            CGM = GM(abs(ii),jj); % CGM=current GM
                 if CGM>TGM
                    bp = [jj ii]; % bp=boundary point
                    Nv = un1';SD=norm(cp-bp); % Nv= normar vector
                    SVMF=1;  %stop vertex moving flag 
                 end
             if norm(bp)==0    
                    if theta1>=-5 && theta1<=5 % inward search start
                            for iii= 1:Dmax
                                 jj=jj+1;
                                CGM = GM(abs(ii),jj); % CGM=current GM
                                if CGM>TGM
                                    bp = [jj ii]; % bp=boundary point
                                    Nv = un1'; SD=norm(cp-bp); % Nv= normar vector                                    
                                    break;
                                end
                            end
                    elseif theta1>=6 && theta1<=84
                           for iii= 1:Dmax
                                 jj=jj+1;ii=ii+1;CGM = GM(abs(ii),jj); 
                                if CGM>TGM
                                    bp = [jj ii]; Nv = un1';SD=norm(cp-bp); 
                                    break;
                                end
                           end
                    elseif theta1>=85 && theta1<=95
                           for iii= 1:Dmax
                                 ii=ii+1;CGM = GM(abs(ii),jj); 
                                if CGM>TGM
                                    bp = [jj ii]; Nv = un1';SD=norm(cp-bp); 
                                    break;
                                end
                           end
                     elseif theta1>=96 && theta1<=174
                           for iii= 1:Dmax
                                 jj=jj-1;ii=ii+1;CGM = GM(abs(ii),jj); 
                                if CGM>TGM
                                    bp = [jj ii]; Nv = un1';SD=norm(cp-bp); 
                                    break;
                                end
                           end  
                     elseif theta1>=175 && theta1<=-175
                           for iii= 1:Dmax
                                 jj=jj-1;CGM = GM(abs(ii),jj); 
                                if CGM>TGM
                                    bp = [jj ii]; Nv = un1'; SD=norm(cp-bp);
                                    break;
                                end
                           end
                     elseif theta1>=-174 && theta1<=-96
                           for iii= 1:Dmax
                                 jj=jj-1; ii=ii-1;CGM = GM(abs(ii),jj);
                                if CGM>TGM
                                    bp = [jj ii];Nv = un1';SD=norm(cp-bp);
                                    break;
                                end
                           end
                     elseif theta1>=-95 && theta1<=-85
                           for iii= 1:Dmax
                                ii=ii-1;CGM = GM(abs(ii),jj);
                                if CGM>TGM
                                    bp = [jj ii];Nv = un1';SD=norm(cp-bp);
                                    break;
                                end
                           end
                     elseif theta1>=-84 && theta1<=-6
                           for iii= 1:Dmax
                                ii=ii-1;jj=jj+1;CGM = GM(abs(ii),jj);
                                if CGM>TGM
                                    bp = [jj ii];Nv = un1';SD=norm(cp-bp);
                                    break;
                                end
                           end
                           %%%%%% outward Search start
                      elseif theta2>=-5 && theta2<=5
                            for iii= 1:Dmax
                                 jj=jj+1;
                                CGM = GM(abs(ii),jj); % CGM=current GM
                                if CGM>TGM
                                    bp = [jj ii]; % bp=boundary point
                                    Nv = run1'; SD=norm(cp-bp); % Nv= normar vector
                                    break;
                                end
                            end
                    elseif theta2>=6 && theta2<=84
                           for iii= 1:Dmax
                                 jj=jj+1;ii=ii+1;CGM = GM(abs(ii),jj); 
                                if CGM>TGM
                                    bp = [jj ii]; Nv = run1';SD=norm(cp-bp); 
                                    break;
                                end
                           end
                    elseif theta2>=85 && theta2<=95
                           for iii= 1:Dmax
                                 ii=ii+1;CGM = GM(abs(ii),jj); 
                                if CGM>TGM
                                    bp = [jj ii]; Nv = run1';SD=norm(cp-bp); 
                                    break;
                                end
                           end
                     elseif theta2>=96 && theta2<=174
                           for iii= 1:Dmax
                                 jj=jj-1;ii=ii+1;CGM = GM(abs(ii),jj); 
                                if CGM>TGM
                                    bp = [jj ii]; Nv = run1';SD=norm(cp-bp); 
                                    break;
                                end
                           end  
                     elseif theta2>=175 && theta2<=-175
                           for iii= 1:Dmax
                                 jj=jj-1;CGM = GM(abs(ii),jj); 
                                if CGM>TGM
                                    bp = [jj ii]; Nv = run1';SD=norm(cp-bp); 
                                    break;
                                end
                           end
                     elseif theta2>=-174 && theta2<=-96
                           for iii= 1:Dmax
                                 jj=jj-1; ii=ii-1;CGM = GM(abs(ii),jj);
                                if CGM>TGM
                                    bp = [jj ii];Nv = run1';SD=norm(cp-bp);
                                    break;
                                end
                           end
                     elseif theta2>=-95 && theta2<=-85
                           for iii= 1:Dmax
                                ii=ii-1;CGM = GM(abs(ii),jj);
                                if CGM>TGM
                                    bp = [jj ii];Nv = run1';SD=norm(cp-bp);
                                    break;
                                end
                           end
                     elseif theta2>=-84 && theta2<=-6
                           for iii= 1:Dmax
                                ii=ii-1;jj=jj+1;CGM = GM(abs(ii),jj);
                                if CGM>TGM
                                    bp = [jj ii];Nv = run1'; SD=norm(cp-bp);
                                    break;
                                end
                           end
                     else
                        SD = 0;
                        Nv = un1';
                    end                                                           
             end  
              %%%%%%%%%%%%%%%%
              
             if SVMF==1
                XX(1,j)=P1(1,1);
                YY(1,j)=P1(1,2);
                VVV(j,1) = 0; % needs to save v1
                VVV(j,2) = 0;
                bp = [0 0];
                SVMF=0;
             else
                wd = 8.9;
                Fd = wd*(SD/Dmax)*Nv; % dynamic distance force
                win = .6; % internal force waiting
                wdamp = -.8;
                wex = .6;
                %F1 = (Fin1*win)+(wex*Fex1);
                F1 = (Fin1*win)+(wex*Fex1)+(wdamp*v1)+Fd;
                %F1 = (wex*Fex1)+(wdamp*v1);
                a1 = (1/m1)*F1;
                v1 = v1+a1;
                VVV(j,1) = v1(1,1); % needs to save v1
                VVV(j,2) = v1(2,1);
                P1 = P1+v1';                      
                XX(1,j)=P1(1,1);
                YY(1,j)=P1(1,2);
                bp = [0 0];
                
             end            
            
    end
X = XX; Y = YY;
% figure(1)
% imagesc(I), colormap(gray(256));
% hold on
% plot(X+0.5,-(Y)+0.5,'b')
% plot(XXX,YYY,'r')
% pause
end
X=X+0.5;
Y=Y+0.5;
Y = -(Y);

figure
imagesc(I), colormap(gray(256));
hold on
plot(X,Y,'b')
plot(XXX,YYY,'r')

figure
imagesc(GM), colormap(gray(256));
hold on
plot(X,Y,'b')
plot(XXX,YYY,'r')




 

    


