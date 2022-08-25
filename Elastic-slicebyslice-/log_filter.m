function Rf = log_filter(R,M)

global filter_params lev

sz = filter_params(lev,1);
sd = filter_params(lev,2);

h = fspecial('gaussian',sz,sd);
     
Rf = padfilter2(h,R);
Rf = Rf.*M;

