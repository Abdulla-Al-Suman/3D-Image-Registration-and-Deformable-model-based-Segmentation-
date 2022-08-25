function [Ny,Nx]= removeoverlap(overlapj,X,Y,L,CNT,NearestMR)

MSRO = 5; %MSRO= Maximum Search remove overlap

iioverlap = abs(Y(overlapj))+0.5; jjoverlap = X(overlapj)+0.5;%jjoverlap= current overlaped point's column

% P1 = [X(overlapj),Y(overlapj)];  % position vectors
%                                        if overlapj==L
%                                         P2 = [X(1),Y(1)];   % Initial nth, 1st, 2nd vertices
%                                        else
%                                         P2 = [X(overlapj+1),Y(overlapj+1)];    
%                                        end
% 
%                                        if overlapj==1
%                                          Pn = [X(L),Y(L)];
%                                       else
%                                         Pn = [X(overlapj-1),Y(overlapj-1)]; 
%                                        end
%                                       d1 = P2-P1;dn=P1-Pn;ud1 = d1/(norm(d1)); udn = dn/(norm(dn));
%                                       c1 = ud1-udn; un1 = c1/(norm(c1));run1 = -(un1);
%                                       theta1 = atan2d(un1(2),un1(1)); % inward searching
%                                       if theta1<0
%                                          theta1 = 360 + theta1;
%                                       end
%                                       theta2 = atan2d(run1(2),run1(1)); % outward searching
%                                        if theta2<0
%                                          theta2 = 360 + theta2;
%                                       end
                                      
                                      CentroidD =CNT-[X(overlapj),Y(overlapj)];UCentroidD = CentroidD/(norm(CentroidD));
                                      theta3 = atan2d(UCentroidD(2),UCentroidD(1));
                                       if theta3<0
                                         theta3 = 360 + theta3;
                                      end
%                                       Dtheta = abs(theta1-theta3); % Dtheta=Difference theta
%                                       if Dtheta<=60
%                                           Stheta = theta1; % Stheta=Search theta
%                                       else
%                                           Stheta = theta2;
%                                       end
                                      %Stheta

                                            if (theta3>=355 && theta3<360) || (theta3>=0 && theta3<5)% along curvature or normal search start
                                                    for io= 1:MSRO
                                                         jjoverlap=jjoverlap+1;
                                                        CGM = NearestMR(round(iioverlap),round(jjoverlap));
                                                        if (CGM==0) || (io== MSRO)  
                                                            Ny = -(iioverlap-0.5);Nx = jjoverlap-0.5;                                    
                                                            break; 
                                                        end
                                                        
                                                    end
                                            elseif theta3>=5 && theta3<85
                                                   for io= 1:MSRO
                                                         jjoverlap=jjoverlap+1;iioverlap=iioverlap-1; 
                                                        CGM = NearestMR(round(iioverlap),round(jjoverlap));
                                                        if (CGM==0) || (io== MSRO) 
                                                            Ny = -(iioverlap-0.5);Nx = jjoverlap-0.5;                                 
                                                            break;
                                                        end
                                                       
                                                   end
                                            elseif theta3>=85 && theta3<95
                                                   for io= 1:MSRO
                                                         iioverlap=iioverlap-1;                               
                                                         CGM = NearestMR(round(iioverlap),round(jjoverlap));
                                                        if (CGM==0) || (io== MSRO) 
                                                            Ny = -(iioverlap-0.5);Nx = jjoverlap-0.5;                                  
                                                            break;
                                                        end
                                                         
                                                   end
                                             elseif theta3>=95 && theta3<175
                                                   for io= 1:MSRO
                                                         jjoverlap=jjoverlap-1;iioverlap=iioverlap-1;                                       
                                                         CGM = NearestMR(round(iioverlap),round(jjoverlap));
                                                        if (CGM==0) || (io== MSRO) 
                                                            Ny = -(iioverlap-0.5);Nx = jjoverlap-0.5;                                  
                                                            break; 
                                                        end
                                                        
                                                   end  
                                             elseif theta3>=175 && theta3<185
                                                   for io= 1:MSRO
                                                         jjoverlap=jjoverlap-1;                                                        
                                                         CGM = NearestMR(round(iioverlap),round(jjoverlap));
                                                        if (CGM==0) || (io== MSRO) 
                                                             Ny = -(iioverlap-0.5);Nx = jjoverlap-0.5;                                   
                                                            break;  
                                                        end
                                                        
                                                   end
                                             elseif theta3>=185 && theta3<265
                                                   for io= 1:MSRO
                                                         jjoverlap=jjoverlap-1; iioverlap=iioverlap+1;                                                    
                                                         CGM = NearestMR(round(iioverlap),round(jjoverlap));
                                                        if (CGM==0) || (io== MSRO)  
                                                           Ny = -(iioverlap-0.5);Nx = jjoverlap-0.5;                                   
                                                            break;  
                                                        end
                                                        
                                                   end
                                             elseif theta3>=265 && theta3<275
                                                   for io= 1:MSRO
                                                        iioverlap=iioverlap+1;                                                                                                           
                                                        CGM = NearestMR(round(iioverlap),round(jjoverlap));
                                                        if (CGM==0) || (io== MSRO)  
                                                            Ny = -(iioverlap-0.5);Nx = jjoverlap-0.5;                                   
                                                            break; 
                                                        end
                                                       
                                                   end
                                             elseif theta3>=275 && theta3<355
                                                   for io= 1:MSRO
                                                        iioverlap=iioverlap+1;jjoverlap=jjoverlap+1;                                                   
                                                       CGM = NearestMR(round(iioverlap),round(jjoverlap));
                                                        if (CGM==0) || (io== MSRO)  
                                                             Ny = -(iioverlap-0.5);Nx = jjoverlap-0.5;                                 
                                                            break; 
                                                        end
                                                        
                                                   end
                                            end