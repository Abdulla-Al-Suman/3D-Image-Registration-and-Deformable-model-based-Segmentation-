function R_hat = calc_R_hat(I0,R0)

% find the size of the images
[n,m] = size(I0);

% set the number of histogram bins
L = 64;

% make sure the values in the image are in the range [0-255]
I0 = (I0<0).*0+(I0>255).*255+((I0>=0)&(I0<=255)).*I0;
R0 = (R0<0).*0+(R0>255).*255+((R0>=0)&(R0<=255)).*R0;

% I0 = (I0<0).*0+floor(I0.*(255/max(max(I0))));
% R0= (R0<0).*0+floor(R0.*(255/max(max(R0))));

I0 = I0-min(min(I0));
I0 = I0.*255./max(max(I0));
R0 = R0-min(min(I0));
R0 = R0.*255./max(max(R0));

% quantize the images
I = round(I0.*L/256);
R = round(R0.*L/256);

% calculate the joint histogram
p = zeros(L+1,L+1);

for u = 1:n
    for v = 1:m
        i = I(u,v);
        r = R(u,v);

            p(i+1,r+1) = p(i+1,r+1)+1;
        
    end
end

% calculate the conditional mean for each value in R
for r = 0:L
    cond_mean(r+1) = sum(p(:,r+1).*[0:L]')./(sum(p(:,r+1))+1e-40);
end

% calculate R_hat by replacing each pixel in R with its conditional mean
% from the joint histogram and inverse quantizing
for u = 1:n
     for v = 1:m
         R_hat(u,v) =  cond_mean(R(u,v)+1).*256/L;
     end
end

