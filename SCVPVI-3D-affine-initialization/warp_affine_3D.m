function V = warp_affine_3D(V0)
global m x y z 

M = [m(1) m(2) m(3) m(4); m(5) m(6) m(7)  m(8);m(9)  m(10) m(11) m(12); 0 0 0 1];
    
xn = M(1,1)*x+M(1,2)*y+M(1,3)*z+M(1,4);
yn = M(2,1)*x+M(2,2)*y+M(2,3)*z+M(2,4);
zn = M(3,1)*x+M(3,2)*y+M(3,3)*z+M(3,4);

%V=interp3(x,y,z,V0,xn,yn,zn,method,0); 
V=interp3_mex(x,y,z,V0,xn,yn,zn); 




