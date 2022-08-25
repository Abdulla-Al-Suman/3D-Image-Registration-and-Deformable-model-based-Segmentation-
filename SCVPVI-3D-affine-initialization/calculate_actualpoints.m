function [xp,yp,zp]=calculate_actualpoints(xp1,yp1,zp1)

 N=256;
 M=256;
 O=45;
 
 x_low=-128;
 y_low=-128;
 z_low=-23;
 

 x_high=127;
 y_high=127;
 z_high=21;
 
 
     s_x = double(1-N)/(x_low - x_high);
     s_y = double(1-M)/(y_low - y_high);
     s_z = double(1-O)/(z_low - z_high);

    o_x = double(1)-x_low.*s_x;
   o_y = double(1)-y_low.*s_y;
     o_z = double(1)-z_low.*s_z;


  xp = s_x*xp1+o_x;
  yp = s_y*yp1+o_y;
  zp = s_z*zp1+o_z;