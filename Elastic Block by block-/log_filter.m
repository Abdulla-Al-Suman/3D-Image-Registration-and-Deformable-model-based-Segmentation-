function Rf = log_filter(R,M)

global filter_params filter_exp lev

sz = filter_params(lev,1);
sd = filter_params(lev,2);

h = fspecial('gaussian',sz,sd);
     
Rf = padfilter2(h,R);
% s = sign(Rf);
% Rf = s.*(abs(Rf).^(filter_exp(lev)));
% Rf = Rf;
Rf = Rf.*M;
% max_Rf = 3*std(Rf(:))
% mu = mean(Rf(:));
% Rf = (Rf-mu).*(104./max_Rf)+mu;

