function V = warp_Elastic_3D(V0,M,s,method)
global x y
[N,MM] = size(V0);
V = zeros(N,MM);

 %[x,y] = meshgrid((0:MM-1)-MM/2,(0:N-1)-N/2);
%[x,y] = meshgrid(0:MM-1,0:N-1);

dxdm = zeros(N,MM,s^2);
k = 1;
for u = 0:s-1
    for v = 0:s-1
        dxdm(:,:,k) = cos(((2*x+1)*pi*u)/(2*MM)).*cos(((2*y+1)*pi*v)/(2*N));
        k = k+1;
    end
end

 xn = zeros(N,MM);
 yn = zeros(N,MM);
%  xn = x;
%  yn = y;
for k = 1:(s^2)
    xn = xn+M(k)*squeeze(dxdm(:,:,k));
    yn = yn+M(k+s^2)*squeeze(dxdm(:,:,k));
end
fxn=xn+x;
fyn=yn+y;



%  V = interp3(x,y,z,V0,xn,yn,zn,'linear',55);
  V = interp2(x,y,V0,fxn,fyn,method);

% figure(2), imagesc([V(:,:,5)  V0(:,:,5)]) 