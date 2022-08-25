function [X,Y] = deformation(TGM,Dmax,X,Y,XX,YY,VVV,k1,dEimdx,dEimdy,GM,wd,win,wex,wdamp,m1,L,bp,SVMF,SD,Nv,NearestMR)

                    for j=1:L % vertex selection loop
                           %j
                            P1 = [X(j),Y(j)];  % position vectors
                                if j==L
                                    P2 = [X(1),Y(1)];   % Initial nth, 1st, 2nd vertices
                                else
                                    P2 = [X(j+1),Y(j+1)];    
                                end

                                if j==1
                                    Pn = [X(L),Y(L)];
                                else
                                    Pn = [X(j-1),Y(j-1)]; 
                                end
                            v1 = [VVV(j,1);VVV(j,2)]; 

                            % calculating Internal forces
                            d1 = P2-P1;dn=P1-Pn;
                            % Resampling 
                    %         ldes = 4; lmin = 0.5*ldes; lmax = 1.5*ldes;
                    %         
                    %         if norm(d1)<lmin
                    %             
                    %         elseif norm(d1)>lmin 
                    %             
                    %         end
                    %         
                            ud1 = d1/(norm(d1)); udn = dn/(norm(dn));
                            ut1 = (ud1+udn)/(norm(ud1+udn)); %unit tangent
                            ur1 = [0 1;-1 0]*ut1'; % column vector ; unit normal
                            c1 = ud1-udn; un1 = c1/(norm(c1)); % un1=unit normal; local curvature
                            run1 = -(un1);% rotated unit normal
                            %ucd1 = (CNT-P1)/(norm(CNT-P1)); rucd1 = -(ucd1);%rucd1= riverse unit centroid direction
                            %k1 = [0 0 0 0 0 -.0005 -.0005 .01 -.0005 -.0005 0 0 0 0 0];
                            %dp = [c1(1)*ur1(1) c1(2)*ur1(2)];
                            Finl1 = conv((dot(c1,ur1')),k1);
                            %Finl1 = conv(dp,k1);
                            %Fin1 = (norm(Finl1))*ur1; % column vector
                            Fin1 = Finl1(8)*ur1; % column vector
                            %Fin1 = (dot(c1,ur1))*ur1;
                            % calculating external force on a vertex
                    %         fx = interp2(xx,yy,dEimdx,X(1),Y(1));
                    %         fy = interp2(xx,yy,dEimdy,X(1),Y(1));
                            fx = interp2(dEimdx,P1(1)+0.5,abs(P1(2))+0.5);
                            fy = interp2(dEimdy,P1(1)+0.5,abs(P1(2))+0.5);
                            Fim1 = [fx,fy];
                            %  z1 = fx+(fy*i);
                            %  r1 = abs(z1);
                            %  theta1 = angle(z1);
                            Fex1 = (dot(Fim1',ur1))*ur1; % column vector

                            %calculating Dynamic distance force
                            theta1 = atan2d(un1(2),un1(1)); % inward searching
                            if theta1<0
                              theta1 = 360 + theta1;
                            end
                            theta2 = atan2d(run1(2),run1(1)); % outward searching
                            if theta2<0
                             theta2 = 360 + theta2;
                            end
        %                     theta1 = atan2d(ucd1(2),ucd1(1)); % inward searching
        %                     theta2 = atan2d(rucd1(2),rucd1(1)); % outward searching
                            ii = ceil(P1(2)); % row
                            jj = ceil(P1(1));
                            cp = [jj ii]; %cp=contour point

                            CGM = GM(abs(ii),jj); % CGM=current GM
                                 if CGM>TGM %&& CGM<TGM2  % checking current point as boundary point
                                    bp = [jj ii]; % bp=boundary point
                                    Nv = [0;0];SD=norm(cp-bp); % Nv= normar vector
                                    SVMF=1;  %stop vertex moving flag 
                                 end

                             if norm(bp)==0 % Starting search

                                            if (theta1>=355 && theta1<360) || (theta1>=0 && theta1<5) % along curvature or normal search start
                                                    for iii= 1:Dmax
                                                         jj=jj+1;
                                                         if (jj>256 || jj<1) || (NearestMR(abs(ii),jj)==1)% for object near to image boundary check
                                                             SD = 0; Nv = [0;0];break;
                                                         end
                                                         
                                                        CGM = GM(abs(ii),jj);
                                                        if CGM>TGM %&& CGM<TGM2
                                                            bp = [jj ii]; % bp=boundary point
                                                            Nv = un1'; SD=norm(cp-bp); % Nv= normar vector                                    
                                                            break;
                                                        end
                                                    end
                                            elseif theta1>=5 && theta1<85
                                                   for iii= 1:Dmax
                                                         jj=jj+1;ii=ii+1;
                                                             if (jj>256 || jj<1) || (NearestMR(abs(ii),jj)==1)
                                                             SD = 0; Nv = [0;0];break;
                                                             elseif abs(ii)>256 || abs(ii)<1
                                                                SD = 0; Nv = [0;0]; break;
                                                              end
                                                         CGM = GM(abs(ii),jj);
                                                        if CGM>TGM %&& CGM<TGM2
                                                            bp = [jj ii]; Nv = un1';SD=norm(cp-bp); 
                                                            break;
                                                        end
                                                   end
                                            elseif theta1>=85 && theta1<95
                                                   for iii= 1:Dmax
                                                         ii=ii+1;
                                                             if (abs(ii)>256 || abs(ii)<1) || (NearestMR(abs(ii),jj)==1)
                                                                SD = 0; Nv = [0;0]; break;
                                                              end
                                                         CGM = GM(abs(ii),jj);
                                                        if CGM>TGM %&& CGM<TGM2
                                                            bp = [jj ii]; Nv = un1';SD=norm(cp-bp); 
                                                            break;
                                                        end
                                                   end
                                             elseif theta1>=95 && theta1<175
                                                   for iii= 1:Dmax
                                                         jj=jj-1;ii=ii+1;
                                                           if (jj>256 || jj<1) || (NearestMR(abs(ii),jj)==1)
                                                             SD = 0; Nv = [0;0];break;
                                                          elseif abs(ii)>256 || abs(ii)<1
                                                                 SD = 0; Nv = [0;0];break;
                                                           end
                                                         CGM = GM(abs(ii),jj);
                                                        if CGM>TGM %&& CGM<TGM2
                                                            bp = [jj ii]; Nv = un1';SD=norm(cp-bp); 
                                                            break;
                                                        end
                                                   end  
                                             elseif theta1>=175 && theta1<185
                                                   for iii= 1:Dmax
                                                         jj=jj-1;
                                                           if (jj>256 || jj<1) || (NearestMR(abs(ii),jj)==1)
                                                            SD = 0; Nv = [0;0]; break;                                                  
                                                           end
                                                         CGM = GM(abs(ii),jj);
                                                        if CGM>TGM %&& CGM<TGM2
                                                            bp = [jj ii]; Nv = un1'; SD=norm(cp-bp);
                                                            break;
                                                        end
                                                   end
                                             elseif theta1>=185 && theta1<265
                                                   for iii= 1:Dmax
                                                         jj=jj-1; ii=ii-1;
                                                           if (jj>256 || jj<1) || (NearestMR(abs(ii),jj)==1)
                                                             SD = 0; Nv = [0;0];break;
                                                           elseif abs(ii)>256 || abs(ii)<1
                                                                SD = 0; Nv = [0;0]; break;
                                                           end
                                                         CGM = GM(abs(ii),jj);
                                                        if CGM>TGM %&& CGM<TGM2
                                                            bp = [jj ii];Nv = un1';SD=norm(cp-bp);
                                                            break;
                                                        end
                                                   end
                                             elseif theta1>=265 && theta1<275
                                                   for iii= 1:Dmax
                                                        ii=ii-1;                                                  
                                                             if (abs(ii)>256 || abs(ii)<1) || (NearestMR(abs(ii),jj)==1)
                                                                SD = 0; Nv = [0;0]; break;
                                                             end
                                                        CGM = GM(abs(ii),jj);
                                                        if CGM>TGM %&& CGM<TGM2
                                                            bp = [jj ii];Nv = un1';SD=norm(cp-bp);
                                                            break;
                                                        end
                                                   end
                                             elseif theta1>=275 && theta1<355
                                                   for iii= 1:Dmax
                                                        ii=ii-1;jj=jj+1;
                                                          if (jj>256 || jj<1) || (NearestMR(abs(ii),jj)==1)
                                                             SD = 0; Nv = [0;0];break;
                                                          elseif abs(ii)>256 || abs(ii)<1
                                                                SD = 0; Nv = [0;0]; break;
                                                          end
                                                        CGM = GM(abs(ii),jj);
                                                        if CGM>TGM %&& CGM<TGM2
                                                            bp = [jj ii];Nv = un1';SD=norm(cp-bp);
                                                            break;
                                                        end
                                                   end
                                                   %%%%%% along reverse curvature or normal Search start
                                              elseif (theta2>=355 && theta2<360) || (theta2>=0 && theta2<5)
                                                    for iii= 1:Dmax
                                                         jj=jj+1;
                                                           if (jj>256 || jj<1) || (NearestMR(abs(ii),jj)==1)
                                                            SD = 0; Nv = [0;0]; break;                                                    
                                                           end
                                                        CGM = GM(abs(ii),jj); 
                                                        if CGM>TGM %&& CGM<TGM2
                                                            bp = [jj ii]; % bp=boundary point
                                                            Nv = run1'; SD=norm(cp-bp); % Nv= normar vector
                                                            break;
                                                        end
                                                    end
                                            elseif theta2>=5 && theta2<85
                                                   for iii= 1:Dmax
                                                         jj=jj+1;ii=ii+1;
                                                           if (jj>256 || jj<1) || (NearestMR(abs(ii),jj)==1)
                                                             SD = 0; Nv = [0;0];break;
                                                           elseif abs(ii)>256 || abs(ii)<1
                                                            SD = 0; Nv = [0;0]; break;
                                                           end
                                                         CGM = GM(abs(ii),jj);
                                                        if CGM>TGM %&& CGM<TGM2
                                                            bp = [jj ii]; Nv = run1';SD=norm(cp-bp); 
                                                            break;
                                                        end
                                                   end
                                            elseif theta2>=85 && theta2<95
                                                   for iii= 1:Dmax
                                                         ii=ii+1;                                                  
                                                           if (abs(ii)>256 || abs(ii)<1) || (NearestMR(abs(ii),jj)==1)
                                                             SD = 0; Nv = [0;0];break;
                                                           end
                                                         CGM = GM(abs(ii),jj);
                                                        if CGM>TGM %&& CGM<TGM2
                                                            bp = [jj ii]; Nv = run1';SD=norm(cp-bp); 
                                                            break;
                                                        end
                                                   end
                                             elseif theta2>=95 && theta2<175
                                                   for iii= 1:Dmax
                                                         jj=jj-1;ii=ii+1;
                                                           if (jj>256 || jj<1) || (NearestMR(abs(ii),jj)==1)
                                                             SD = 0; Nv = [0;0];break;
                                                           elseif abs(ii)>256 || abs(ii)<1
                                                             SD = 0; Nv = [0;0];break;
                                                           end
                                                         CGM = GM(abs(ii),jj);
                                                        if CGM>TGM %&& CGM<TGM2
                                                            bp = [jj ii]; Nv = run1';SD=norm(cp-bp); 
                                                            break;
                                                        end
                                                   end  
                                             elseif theta2>=175 && theta2<185
                                                   for iii= 1:Dmax
                                                         jj=jj-1;
                                                           if (jj>256 || jj<1) || (NearestMR(abs(ii),jj)==1)
                                                             SD = 0; Nv = [0;0];break;                                                     
                                                           end
                                                         CGM = GM(abs(ii),jj); 
                                                        if CGM>TGM %&& CGM<TGM2
                                                            bp = [jj ii]; Nv = run1';SD=norm(cp-bp); 
                                                            break;
                                                        end
                                                   end
                                             elseif theta2>=185 && theta2<265
                                                   for iii= 1:Dmax
                                                         jj=jj-1; ii=ii-1;
                                                           if (jj>256 || jj<1) || (NearestMR(abs(ii),jj)==1)
                                                             SD = 0; Nv = [0;0];break;
                                                           elseif abs(ii)>256 || abs(ii)<1
                                                             SD = 0; Nv = [0;0];break;
                                                           end
                                                         CGM = GM(abs(ii),jj);
                                                        if CGM>TGM %&& CGM<TGM2
                                                            bp = [jj ii];Nv = run1';SD=norm(cp-bp);
                                                            break;
                                                        end
                                                   end
                                             elseif theta2>=265 && theta2<275
                                                   for iii= 1:Dmax
                                                        ii=ii-1;
                                                             if (abs(ii)>256 || abs(ii)<1) || (NearestMR(abs(ii),jj)==1)
                                                                 SD = 0; Nv = [0;0];
                                                                 break;
                                                              end
                                                        CGM = GM(abs(ii),jj);
                                                        if CGM>TGM %&& CGM<TGM2
                                                            bp = [jj ii];Nv = run1';SD=norm(cp-bp);
                                                            break;
                                                        end
                                                   end
                                             elseif theta2>=275 && theta2<355
                                                   for iii= 1:Dmax
                                                        ii=ii-1;jj=jj+1;
                                                          if (jj>256 || jj<1) || (NearestMR(abs(ii),jj)==1)
                                                              SD = 0; Nv = [0;0];
                                                             break;
                                                          elseif abs(ii)>256 || abs(ii)<1
                                                              SD = 0; Nv = [0;0];
                                                              break;
                                                          end
                                                        CGM = GM(abs(ii),jj);
                                                        if CGM>TGM %&& CGM<TGM2
                                                            bp = [jj ii];Nv = run1'; SD=norm(cp-bp);
                                                            break;
                                                        end
                                                   end
                                             else
                                                SD = 0; Nv = [0;0];
                                            end 

                             end % Distance sersing loop 
                              %%%%%%%%%%%%%%%%

                             if SVMF==1
                                XX(1,j)=P1(1,1);YY(1,j)=P1(1,2);
                                VVV(j,1) = 0;VVV(j,2) = 0;bp = [0 0]; SVMF=0;
                             else                      
                                Fd = wd*(SD/Dmax)*Nv; % dynamic distance force                      
                                %F1 = (Fin1*win)+(wex*Fex1)+Fd;
                                F1 = (Fin1*win)+(wex*Fex1)+(wdamp*v1)+Fd;
                                %F1 = (wex*Fex1)+(wdamp*v1);
                                a1 = (1/m1)*F1; v1 = v1+a1;
                                VVV(j,1) = v1(1,1); % needs to save v1
                                VVV(j,2) = v1(2,1);
                                P1 = P1+v1';                      
                                XX(1,j)=P1(1,1);YY(1,j)=P1(1,2);bp = [0 0];%SD = 0; Nv = [0;0];
                             end            
                    end %  Vertex selection loop
                X = XX; Y = YY;
 
