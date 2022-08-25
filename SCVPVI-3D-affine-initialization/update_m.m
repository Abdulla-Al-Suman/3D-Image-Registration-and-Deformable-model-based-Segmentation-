
function update_m(dIdm,I,R)

global m
%[R_hat] = calc_R_hat_pvi_mex(I,R,m);

[R_hat] = calc_R_hat_pvi_new(I,R);

E = I-R_hat;

out1 = 0.0;
out2 = 0.0;

for i = 1:12
    for j = 1:12

        %H(i,j) = sum(sum(sum(dIdm(:,:,:,i).*dIdm(:,:,:,j))));        
        H(i,j)=mul_hes_mex(dIdm(:,:,:,i),dIdm(:,:,:,j),out1);

    end
        %b(i) = sum(sum(sum(squeeze(E.*squeeze(dIdm(:,:,:,i))))));
        b(i) = mul_hes_mex(E,dIdm(:,:,:,i),out2);

end

dm = -inv(H)*b';
%ddm= dm*1.5;
m = m+dm';
    
