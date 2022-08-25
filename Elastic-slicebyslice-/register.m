function [m,I] = register(m,I0,R0,s)

global sm lev num_iter
[h,w] = size(I0);
M0 = ones(h,w);

I = warp_Elastic_3D(I0,m,s,'linear');
M = warp_Elastic_3D(M0,m,s,'linear')==1;

% apply LoG filter to fluoro (R0) and CT projection (I0) 
If = log_filter(I,M);
Rf = log_filter(R0,M);

for iter = 1:num_iter(lev)
    
    % calculate gradients with respect to motion parameters
    dIdm = calc_dIdm(If,s);    

    % find updated motion parameters using selected similarity measure
    if sm == 0
        m = update_m(dIdm,If,Rf,m,s);
    else
        m = update_m_MI(dIdm,If,Rf,m);
    end

    % place I0 at updated position
I = warp_Elastic_3D(I0,m,s,'linear');
M = warp_Elastic_3D(M0,m,s,'linear')==1;

    % apply LoG filter updated CT projection (I0) 
    If = log_filter(I,M);
    Rf = log_filter(R0,M);
end



function W = edge_mask(bdr)

W = [zeros(1,bdr) ones(1,128-2*bdr) zeros(1,bdr)]'*[zeros(1,bdr) ones(1,128-2*bdr) zeros(1,bdr)];



function dIdm = calc_dIdm(I,s)
global x y 
 
[nn,mm] = size(I);
 
[dIdx,dIdy] = gradient(I);
 
dxdm = zeros(nn,mm,s^2);
k = 1;
for u = 0:s-1
    for v = 0:s-1
        dxdm(:,:,k) = cos(((2*x+1)*pi*u)/(2*mm)).*cos(((2*y+1)*pi*v)/(2*nn));
        k = k+1;
    end
end

dIdm = zeros(nn,mm,2*s^2);
k = 1;
    for u = 0:s-1
        for v = 0:s-1

            dIdm(:,:,k) = dIdx.*squeeze(dxdm(:,:,k));%+ dIdy.*squeeze(dxdm(:,:,k));
            dIdm(:,:,k+s^2) = dIdy.*squeeze(dxdm(:,:,k));

            k = k+1;
        end
    end


function V = update_V(m,V0)

global xs ys cog 

M_persp_z = makeM(0,0,0,0,0,0,0,0,0,0,0,0.0014,128-xs,128-ys,64);
M_rigid = makeM(m(1),m(2),m(3),0,0,0,m(4),m(5),m(6),0,0,0,cog(2),cog(1),cog(3));

S = makeresampler('cubic','fill');
T = maketform('projective',(M_persp_z*M_rigid)');

V = tformarray(V0,T,S,[1 2 3],[1 2 3],[128 128 128],[],[]);






