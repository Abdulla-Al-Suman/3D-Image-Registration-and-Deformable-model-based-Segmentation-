function [R0,FI] = test_range(m,s)

global Im lev V VV MFEF 

I0 = zeros(256,256,45);
I0(:,:,:) = V;
I0 = floor(I0.*(255/max(max(max(I0)))));

% Io have to divided for local registration
[M,N,OO] = size(I0);

I = zeros(128,128,180);
pp = 1;
for i = 1:128:(M-127)
  for j = 1:128:(N-127)
      for k = 1:OO
    tmp = I0(i:(i+127), j:(j+127),k);
    I(:,:,pp) = tmp;
    pp = pp+1;
    
      end
  end
end


%load('P2_16.mat')
R0 = zeros(256,256,45);
R0(:,:,:) = VV;
R0= floor(R0.*(255/max(max(max(R0)))));
% R = 256-R;


% Ro have to divided for local registration

R = zeros(128,128,180);    
pp = 1;
for i = 1:128:(M-127)
  for j = 1:128:(N-127)
      for k = 1:OO
    tmp = R0(i:(i+127), j:(j+127),k);
    R(:,:,pp) = tmp;
    pp = pp+1;
      end
  end
end 


FI = zeros(M,N,OO);
% apply registration for each block
    
  i = 1;
    for rs = 1:2
        for cs = 1:2
            for hs = 1:OO

                for lev = 1:3
                    [m,RI] = register(m,squeeze(I(:,:,i)),squeeze(R(:,:,i)),s);
                end
            
                FI((rs-1)*128+1:rs*128,(cs-1)*128+1:cs*128,hs) = RI;
                MFEF(i,:) = m;
                i = i+1;
                m = Im;
            end

        end
    end

