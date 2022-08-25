function m = update_m(dIdm,If,Rf,m,s)
global PDF
Rf_hat = calc_R_hat(If,Rf);

E = If-Rf_hat;

for i = 1:2*s^2
    for j = 1:2*s^2
        H(i,j) = sum(sum(squeeze(dIdm(:,:,i)).*squeeze(dIdm(:,:,j))));
    end
    b(i) = sum(sum(squeeze(E.*squeeze(dIdm(:,:,i)))));
end

dm = -inv(H)*b';
ddm = dm/PDF;

m = m+ddm';
    
