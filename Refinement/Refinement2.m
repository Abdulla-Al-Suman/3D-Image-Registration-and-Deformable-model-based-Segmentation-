clc
clear all 
close all

load('p2 muscle data.mat')   % muscle_data=cell array
load('FR881832-80808060m10304050.mat')
I = zeros(256,256);
I(:,:) = V(:,:,3);
% [D1,D2] = size(I);
% [xx,yy] = meshgrid([1:D2]-120,[1:D1]-123);
%a = min(I)

% Contour x y values retriving
X = muscle_data{5}{1}{3}.x;  % 5th muscle 1= left side, 1=slice x data
Y = muscle_data{5}{1}{3}.y;  % Y= y coordinates values for all vertices
XXX = X; YYY = Y; % saving initial positions
X=X-0.5; Y=Y-0.5;
Y = -(Y);
L = length(X);



% XX = zeros(L,2); % XX= initial contour matrix
% for i = 1:L
%     XX(i,1) = X(i);
%     XX(i,2) = Y(i);
% end

% Calculating potential energy function/distribution
G=fspecial('gauss',[20 20],5);
CI(:,:)=filter2(G,I);
[dCIdx,dCIdy] = gradient(CI);

dx = zeros(256,256);
dy = zeros(256,256);

for i=1:256
    for j=1:256
        dx(i,j) = dCIdx(i,j)^2;
        dy(i,j) = dCIdy(i,j)^2;
    end
end
Eim = -(dx+dy); % - means valley
%Eim = 0.5*CI; % - means valley
%imagesc(Eim), colormap(gray(256));

[dEimdx,dEimdy] = gradient(Eim);

dEimdx = -1*dEimdx; % external force field
dEimdy = -1*dEimdy; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


XX = zeros(1,L); % deformed contour x values for all vertices
YY = zeros(1,L);
VVV = zeros(L,2); % initial velocities
% AA = zeros(1,L); % for acceleration saving
for i=1:80
    

    for j=1:L

        P1 = [X(j),Y(j)];  % position vectors

        if j==L
            P2 = [X(1),Y(1)];   % Initial nth 1st 2nd vertices
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
        %norm(ud1)

        ut1 = (ud1+udn)/(norm(ud1+udn));

        ur1 = [0 1;-1 0]*ut1'; % column vector

        c1 = ud1-udn;

        k1 = [0 0 0 0 0 0 -.5 2 -.5 0 0 0 0 0 0];
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
        %%%%%%%%%%%%%%%%


        win = .05; % internal force waiting
        wdamp = -.7;
        wex = .6;
        %F1 = (Fin1*win)+(wex*Fex1);
        F1 = (Fin1*win)+(wex*Fex1)+(wdamp*v1);
        m1 = 1;
        a1 = (1/m1)*F1;
        P1 = P1+v1';
%         AA(j) = norm(a1);
        v1 = v1+a1;
        VVV(j,1) = v1(1,1); % needs to save v1
        VVV(j,2) = v1(2,1);
        XX(1,j)=P1(1,1);
        YY(1,j)=P1(1,2);
    end
 X = XX; Y = YY;
 
 % v,a=0 check
%  if norm(AA)==0
%     break
%  end

     
end 
X=X+0.5; Y=Y+0.5;
Y = -(Y);
figure
imagesc(I), colormap(gray(256));
hold on
plot(X,Y,'b')
plot(XXX,YYY,'r')

figure
imagesc(Eim), colormap(gray(256));
hold on
plot(X,Y,'b')




 

    


