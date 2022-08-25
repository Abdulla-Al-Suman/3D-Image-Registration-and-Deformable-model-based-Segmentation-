clc
clear all 
close all
load('R-muscleDataP2toP1.mat') % contains autosegmentation results
load('Patient 1 muscle data.mat') % contains ground truth
load('T-muscleDataP2toP1.mat') % for transfer DSC calculation
DSCtableR = zeros(12,12); % for containing DSC value
DSCtableT = zeros(12,12);
DSCtableRS = zeros(12,12);DSCtableTS = zeros(12,12);
musclesmeanR = zeros(1,12); % FOR CONTAINING MEAN VALUES FOR EACH MUSCLES
musclesmeanT = zeros(1,12);
musclesmeanTS = zeros(1,12);musclesmeanRS = zeros(1,12);
musclesvarianceTS = zeros(1,12);musclesvarianceRS = zeros(1,12);


muscle = [3 11 16 19 22 24];

% Initial DSC table creation
for MN = 1:6
    for SN = 1:2  % 1=left
        for FN = 14:25
            MX = muscle_data{muscle(MN)}{SN}{FN}.x;MY = muscle_data{muscle(MN)}{SN}{FN}.y;LM = length(MX); % ground truth
            X = R_muscle_data{muscle(MN)}{SN}{FN}.x;Y = R_muscle_data{muscle(MN)}{SN}{FN}.y; L = length(X); % Auto segmented musk
            XT = T_muscle_data{muscle(MN)}{SN}{FN}.x;YT = T_muscle_data{muscle(MN)}{SN}{FN}.y;LT = length(XT);
            if L>=1 && LM>=1
            % complete boundary plot
            CX=zeros(1,L+1);CX(1,1:L)=X;CX(1,L+1)=X(1,1);CY=zeros(1,L+1);CY(1,1:L)=Y;CY(1,L+1)=Y(1,1);
            CXT=zeros(1,LT+1);CXT(1,1:LT)=XT;CXT(1,LT+1)=XT(1,1);CYT=zeros(1,LT+1);CYT(1,1:LT)=YT;CYT(1,LT+1)=YT(1,1);
            CMX=zeros(1,LM+1);CMX(1,1:LM)=MX;CMX(1,LM+1)=MX(1,1);CMY=zeros(1,LM+1);CMY(1,1:LM)=MY;CMY(1,LM+1)=MY(1,1);
            
            A = poly2mask(CX, CY, 256, 256);AT = poly2mask(CXT, CYT, 256, 256);M = poly2mask(CMX, CMY, 256, 256);
            AiM = and(A,M);AiMT = and(AT,M);A = double(A);AT = double(AT);M = double(M);AiM = double(AiM);AiMT = double(AiMT);
            NA= sum(A(:));NAT= sum(AT(:)); NM= sum(M(:));NAM=sum(AiM(:));NAMT=sum(AiMT(:));
            DSCtableR((FN-13),(MN-1)*2+SN) = ((2*NAM)/(NA+NM)); % need to change
            DSCtableT((FN-13),(MN-1)*2+SN) = ((2*NAMT)/(NAT+NM)); 
            end
            
        end
    end
end
% DSC selected table creation
for MN = 1:12
    for FN = 1:12        
        if DSCtableR(FN,MN)>0.75  
         DSCtableRS(FN,MN) = DSCtableR(FN,MN);
        end
        if DSCtableT(FN,MN)>0.75  
         DSCtableTS(FN,MN) = DSCtableT(FN,MN);
        end
    end
end

% Muscle Mean calculation
i=0;SmuscleR=0;j=0;SmuscleT=0;k=0;SmuscleTS=0;l=0;SmuscleRS=0;
 VarianceATS = zeros(12,1); VarianceARS = zeros(12,1);%VarianceATS= variance array transformed selected 

for MN = 1:12
    for FN = 1:12
        if DSCtableR(FN,MN)>0
            i=i+1;SmuscleR= SmuscleR+DSCtableR(FN,MN);
            j=j+1;SmuscleT= SmuscleT+DSCtableT(FN,MN);
        end
        if DSCtableTS(FN,MN)>0.80
            k=k+1;SmuscleTS= SmuscleTS+DSCtableTS(FN,MN);
            VarianceATS(k,1) = DSCtableTS(FN,MN);
        end
        if DSCtableRS(FN,MN)>0.80
            l=l+1;SmuscleRS= SmuscleRS+DSCtableRS(FN,MN);
            VarianceARS(l,1) = DSCtableRS(FN,MN);
        end
    end
    musclesmeanR(MN)=SmuscleR/i;i=0;SmuscleR=0;
    musclesmeanT(MN)=SmuscleT/j;j=0;SmuscleT=0;
    musclesmeanTS(MN)=SmuscleTS/k;k=0;SmuscleTS=0;
    musclesmeanRS(MN)=SmuscleRS/l;l=0;SmuscleRS=0;
    musclesvarianceTS(MN) = var(VarianceATS);VarianceATS=0;
    musclesvarianceRS(MN) = var(VarianceARS);VarianceARS=0;
end
DSCtableT
DSCtableR
DSCtableTS
DSCtableRS
musclesmeanT
musclesmeanR
musclesmeanTS
musclesmeanRS
musclesvarianceTS
musclesvarianceRS


TmeanRS =  zeros(1,8);TmeanTS =  zeros(1,8);% Tmean= total mean
for MN=1:6
    TmeanRS(MN)= musclesmeanRS(MN);
    TmeanTS(MN)= musclesmeanTS(MN);
end
for MN=9:10
    TmeanRS(MN-2)= musclesmeanRS(MN);
    TmeanTS(MN-2)= musclesmeanTS(MN);
end
meanDSCRS = mean(TmeanRS)
meanDSCTS = mean(TmeanTS)
TvarienceRS = var(TmeanRS)
TvarienceTS = var(TmeanTS)

% meanDSCRS = mean(musclesmeanR)
% meanDSCTS = mean(musclesmeanT)
%meanDSCTS = mean(musclesmeanTS)



