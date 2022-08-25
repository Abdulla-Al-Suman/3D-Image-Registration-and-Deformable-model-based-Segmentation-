function VVV = warp_Elastic_3D(V0,M,s,method)
global x y 
[N,MM] = size(V0);
VVV = zeros(N,MM);

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
for k = 1:(s^2)
    xn = xn+M(k)*squeeze(dxdm(:,:,k));
    yn = yn+M(k+s^2)*squeeze(dxdm(:,:,k));
end
fxn=xn+x;
fyn=yn+y;
VVV = interp2(x,y,V0,fxn,fyn,method,0);
 