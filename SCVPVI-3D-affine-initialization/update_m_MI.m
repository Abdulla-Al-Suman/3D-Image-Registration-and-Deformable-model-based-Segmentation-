function m = update_m_MI(dIdm,I,R,m)


[p,pt,dpdm,dptdm,dSdm] = calc_dSdm_MI(I,R,dIdm);


    
for i = 1:6
    for j = 1:6
        H(i,j) = 1.4427.*(sum(squeeze(dptdm(i,:)).*squeeze(dptdm(j,:)).*(pt~=0)./(pt+1e-40))...
            -sum(sum(squeeze(dpdm(i,:,:)).*squeeze(dpdm(j,:,:)).*(p~=0)./(p+1e-40))));
    end
    b(i) = -sum(sum(squeeze(squeeze(dSdm(i,:,:)))));
end

dm = -inv(H)*b';

m = m+25.*dm';

