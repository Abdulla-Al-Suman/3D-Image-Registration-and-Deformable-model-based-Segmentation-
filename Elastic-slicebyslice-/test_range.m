function [R0,FI] = test_range(m,s)

global lev Im V VV MFEF 

I0 = zeros(256,256,45);
I0(:,:,:) = V;
I0 = floor(I0.*(255/max(max(max(I0)))));
[M,N,OO] = size(I0);

R0 = zeros(256,256,45);
R0(:,:,:) = VV;
R0= floor(R0.*(255/max(max(max(R0)))));

FI = zeros(M,N,OO);
% apply registration for each block
    
% i = 1;
%     for rs = 1:4
        for zs = 1:OO
            for lev = 1:3
                [m,RI] = register(m,squeeze(I0(:,:,zs)),squeeze(R0(:,:,zs)),s);
            end            
                FI(:,:,zs) = RI;
                MFEF(zs,:) = m; % MFEF = m for each frame
                m = Im;
                zs              
        end
%     end

