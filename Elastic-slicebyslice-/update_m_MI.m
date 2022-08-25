function m = update_m_CCRE(dIdm,If,Rf,m)

global lev CCRE_step_size

[p,pt,dpdm,dptdm,dSdm] = calc_dSdm_MI(If,Rf,dIdm);


p1 = [1 2 6];
p2 = [1 2 4 5 6];

a = CCRE_step_size(lev);

    
for i = 1:6
    for j = 1:6
        H(i,j) = 1.4427.*(sum(squeeze(dptdm(i,:)).*squeeze(dptdm(j,:)).*(pt~=0)./(pt+1e-40))...
            -sum(sum(squeeze(dpdm(i,:,:)).*squeeze(dpdm(j,:,:)).*(p~=0)./(p+1e-40))));
    end
    b(i) = -sum(sum(squeeze(squeeze(dSdm(i,:,:)))));
end

dm = -inv(H)*b';

m = m+a.*dm';

