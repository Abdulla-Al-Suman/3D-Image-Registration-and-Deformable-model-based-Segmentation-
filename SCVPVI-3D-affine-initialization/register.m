function [I] = register(I0,R0)

global sm  lev num_iter m 

% place I0 and M0 at initial position

[D1,D2,D3] = size(I0);
M0 = ones(D1,D2,D3);

I = warp_affine_3D(I0);
M = warp_affine_3D(M0)==1;

for iter = 1:num_iter(lev)
    % calculate gradients with respect to motion parameters
    dIdm = calc_dIdm(I);

    % find updated motion parameters using selected similarity measure
    if sm == 0
        update_m(dIdm,I,R0);
    else
        m = update_m_MI(dIdm,I.*M,R0.*M,m);
    end
    % place I0 at updated position
     I = warp_affine_3D(I0);    
     M = warp_affine_3D(M0)==1;
end


function dIdm = calc_dIdm(I)
global x y z

[D1,D2,D3] = size(I);

[dIdx,dIdy,dIdz] = gradient(I);

% Find values for dI'/dm1, dI'/dm2, dI'/dm3 and  dI'/dm4

dIdm(:,:,:,1) = x.*dIdx;
dIdm(:,:,:,2) = y.*dIdx;
dIdm(:,:,:,3) = z.*dIdx;
dIdm(:,:,:,4) = dIdx;

% Find values for dI'/dm5, dI'/dm6, dI'/dm7 and dI'/dm8
dIdm(:,:,:,5) = x.*dIdy;
dIdm(:,:,:,6) = y.*dIdy;
dIdm(:,:,:,7) = z.*dIdy;
dIdm(:,:,:,8) = dIdy;

% Find values for dI'/dm9, dI'/dm10, dI'/dm11 and dI'/dm12
dIdm(:,:,:,9) = x.*dIdz;
dIdm(:,:,:,10) = y.*dIdz;
dIdm(:,:,:,11) = z.*dIdz;
dIdm(:,:,:,12) = dIdz;

