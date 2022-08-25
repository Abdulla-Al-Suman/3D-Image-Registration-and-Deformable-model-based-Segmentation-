function y = padfilter2(b,x)

[n,m] = size(x);
[u,v] = size(b);

padx = zeros(n+2*u,m+2*v);
xf = zeros(n+2*u,m+2*v);
y = zeros(n,m);

padx(u+1:n+u,v+1:m+v) = x;
padx(1:u,:) = imresize_old(padx(u+2,:),[u m+2*v]);
padx(n+u+1:n+2*u,:) = imresize_old(padx(n+u-1,:),[u m+2*v]);
padx(:,1:v) = imresize_old(padx(:,v+2),[n+2*u v]);
padx(:,m+v+1:m+2*v) = imresize_old(padx(:,m+v-1),[n+2*u v]);

xf = filter2(b,padx);

y = xf(u+1:n+u,v+1:m+v);
