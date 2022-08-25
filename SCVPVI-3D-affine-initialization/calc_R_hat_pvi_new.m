function [R_hat] = calc_R_hat_pvi_new(I0,R0)
global m
% find the size of the images
[d1,d2,d3] = size(I0);
% set the number of histogram bins
L =256;

% make sure the values inthe image are in the range [0-255]
I0 = (I0<0).*0+(I0>255).*255+((I0>=0)&(I0<=255)).*I0;
R0 = (R0<0).*0+(R0>255).*255+((R0>=0)&(R0<=255)).*R0;

I0 = I0-min(min(min(I0)));
I0 = I0.*255./max(max(max(I0)));
R0 = R0-min(min(min(I0)));
R0 = R0.*255./max(max(max(R0)));
 
% quantize the images
I = round(I0.*L/256);
R = round(R0.*L/256);


M = [m(1) m(2) m(3) m(4); m(5) m(6) m(7)  m(8);m(9)  m(10) m(11) m(12); 0 0 0 1];
    
[x1,y1,z1] = meshgrid([0:d2-1]-d2/2,[0:d1-1]-d1/2,[0:d3-1]-d3/2);

%[x1,y1,z1] = meshgrid([0:d2-1],[0:d1-1],[0:d3-1]);
x = M(1,1)*x1+M(1,2)*y1+M(1,3)*z1+M(1,4);
y = M(2,1)*x1+M(2,2)*y1+M(2,3)*z1+M(2,4);
z = M(3,1)*x1+M(3,2)*y1+M(3,3)*z1+M(3,4);

% To calculate weight and joint-hist for the images

[y,x,z]=calculate_actualpoints(y,x,z);
              

h=zeros(256+1,256+1);


for i=1:d1
 
      for j=1:d2
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
          for k=1:d3
              
              
              %xp1=y(i,j,k); yp1=x(i,j,k); zp1=z(i,j,k); % Actual points
       
             % [xp,yp,zp]=calculate_actualpoints(xp1,yp1,zp1);
              
             xp=y(i,j,k); yp=x(i,j,k); zp=z(i,j,k);
             
              fx=floor(xp); fy=floor(yp); fz=floor(zp); 
              cx=floor(xp)+1; cy=floor(yp)+1; cz=floor(zp)+1;  % Ceil is equivalent to floor+1
        
     
              %if (xp>=1 && xp<=256) && (yp>=1 && yp<=256) && (zp>=1 && zp<=32)
              if (fx>=1 && cx<=256) && (fy>=1 && cy<=256) && (fz>=1 && cz<=45)
              
               w1=(cx-xp)*(cy-yp)*(cz-zp);
               w2=(cx-xp)*(cy-yp)*(zp-fz);
               w3=(cx-xp)*(yp-fy)*(cz-zp);
               w4=(cx-xp)*(yp-fy)*(zp-fz);
               w5=(xp-fx)*(cy-yp)*(cz-zp);
               w6=(xp-fx)*(cy-yp)*(zp-fz);
               w7=(xp-fx)*(yp-fy)*(cz-zp);                          
               w8=(xp-fx)*(yp-fy)*(zp-fz);
                             
               % pixel intensities for caculating joint-hist
               v=I(i,j,k);
                                                           
               u1=R(fx,fy,fz);
               u2=R(fx,fy,cz);
               u3=R(fx,cy,fz);
               u4=R(fx,cy,cz);
               u5=R(cx,fy,fz);
               u6=R(cx,fy,cz);
               u7=R(cx,cy,fz);
               u8=R(cx,cy,cz);
               

               h(v+1,u1+1)=h(v+1,u1+1)+w1;
               h(v+1,u2+1)=h(v+1,u2+1)+w2;
               h(v+1,u3+1)=h(v+1,u3+1)+w3;
               h(v+1,u4+1)=h(v+1,u4+1)+w4;
               h(v+1,u5+1)=h(v+1,u5+1)+w5;
               h(v+1,u6+1)=h(v+1,u6+1)+w6;
               h(v+1,u7+1)=h(v+1,u7+1)+w7;
               h(v+1,u8+1)=h(v+1,u8+1)+w8;
                   
             
              end              
          
          end
    
      end
end
% calculate the conditional mean for each value in R
cond_mean = zeros(1,L+1);
for r = 0:L
    cond_mean(r+1) = sum(h(:,r+1).*[0:L]')./(sum(h(:,r+1))+1e-40);   
    %cond_mean(r+1)=max(h(:,r+1));
end


% calculate R_hat by replacing each pixel in R with its conditional mean
% from the joint histogram and inverse quantizing

R_hat = zeros(d1,d2,d3);

for u = 1:d1
    for v = 1:d2
        for w = 1:d3
            
            R_hat(u,v,w) =  cond_mean(R(u,v,w)+1);
        end
    end
end