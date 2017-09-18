%Performs simulations for the case "binary; approx. linear; majority rule; binary"
%AND
%Performs simulations for the case "binary; approx. linear; CDT; binary"



%load EMG dataset and set HD parameters
%clear;
load('COMPLETE_cleaned.mat')
holoGN; %preload functions used during simulations
learningFrac = 0.25;
D = 10000;
percision = 1;
cuttingAngle = 0.9;



global encoderMode;
%encoderMode = 'N_feature_perm';
encoderMode = '4N_features';


global codingMode;
%codingMode = 'dense_binary';
%codingMode = 'dense_bipolar';

global trainingMode;
%trainingMode = 'binary_vect';
%trainingMode = 'integer_vect';
global bipolarDistretizationMode;
bipolarDistretizationMode = 'off';

st = 2;
res = 1;
sp = 20;

downSampRate = 250;
overlap = 0;
[TS_COMPLETE_1, L_TS_COMPLETE_1] = downSampling (COMPLETE_1, LABEL_1, downSampRate);
[TS_COMPLETE_2, L_TS_COMPLETE_2] = downSampling (COMPLETE_2, LABEL_2, downSampRate);
[TS_COMPLETE_3, L_TS_COMPLETE_3] = downSampling (COMPLETE_3, LABEL_3, downSampRate);
[TS_COMPLETE_4, L_TS_COMPLETE_4] = downSampling (COMPLETE_4, LABEL_4, downSampRate);
downSampRate = 50;
[TS_COMPLETE_5, L_TS_COMPLETE_5] = downSampling (COMPLETE_5, LABEL_5, downSampRate);

%DK. Changed argument in getTrainData from '-------' to 'inorder'. See comment in
%runHoloGN. Note that now Training data is always the same
%DK. Manual division into training and testing parts. The training part
%includes 1st, 2nd, and 10th gestures for eacch type

%subject 1
g1_tr_2_end=29;
g1_tr_10_start=127;
g2_tr_2_end=174;
g2_tr_10_start=274;
g3_tr_2_end=320;
g3_tr_10_start=427;
g4_tr_2_end=472;
g4_tr_10_start=575;

SAMPL_DATA_1= [TS_COMPLETE_1(1:g1_tr_2_end,:);TS_COMPLETE_1(g1_tr_10_start:g2_tr_2_end,:); TS_COMPLETE_1(g2_tr_10_start:g3_tr_2_end,:); TS_COMPLETE_1(g3_tr_10_start:g4_tr_2_end,:); TS_COMPLETE_1(g4_tr_10_start:end,:);  ];
L_SAMPL_DATA_1= [L_TS_COMPLETE_1(1:g1_tr_2_end,:);L_TS_COMPLETE_1(g1_tr_10_start:g2_tr_2_end,:); L_TS_COMPLETE_1(g2_tr_10_start:g3_tr_2_end,:); L_TS_COMPLETE_1(g3_tr_10_start:g4_tr_2_end,:); L_TS_COMPLETE_1(g4_tr_10_start:end,:);];
%assign new values to TS_COMPLETE and L_TS_COMPLETE
TS_COMPLETE_1=[TS_COMPLETE_1(g1_tr_2_end+1:g1_tr_10_start-1,:); TS_COMPLETE_1(g2_tr_2_end+1:g2_tr_10_start-1,:); TS_COMPLETE_1(g3_tr_2_end+1:g3_tr_10_start-1,:); TS_COMPLETE_1(g4_tr_2_end+1:g4_tr_10_start-1,:); ];
L_TS_COMPLETE_1=[L_TS_COMPLETE_1(g1_tr_2_end+1:g1_tr_10_start-1,:); L_TS_COMPLETE_1(g2_tr_2_end+1:g2_tr_10_start-1,:); L_TS_COMPLETE_1(g3_tr_2_end+1:g3_tr_10_start-1,:); L_TS_COMPLETE_1(g4_tr_2_end+1:g4_tr_10_start-1,:); ];

%subject 2
g1_tr_2_end=44;
g1_tr_10_start=129;
g2_tr_2_end=179;
g2_tr_10_start=264;
g3_tr_2_end=313;
g3_tr_10_start=394;
g4_tr_2_end=448;
g4_tr_10_start=531;

SAMPL_DATA_2= [TS_COMPLETE_2(1:g1_tr_2_end,:);TS_COMPLETE_2(g1_tr_10_start:g2_tr_2_end,:); TS_COMPLETE_2(g2_tr_10_start:g3_tr_2_end,:); TS_COMPLETE_2(g3_tr_10_start:g4_tr_2_end,:); TS_COMPLETE_2(g4_tr_10_start:end,:);  ];
L_SAMPL_DATA_2= [L_TS_COMPLETE_2(1:g1_tr_2_end,:);L_TS_COMPLETE_2(g1_tr_10_start:g2_tr_2_end,:); L_TS_COMPLETE_2(g2_tr_10_start:g3_tr_2_end,:); L_TS_COMPLETE_2(g3_tr_10_start:g4_tr_2_end,:); L_TS_COMPLETE_2(g4_tr_10_start:end,:);];
%assign new values to TS_COMPLETE and L_TS_COMPLETE
TS_COMPLETE_2=[TS_COMPLETE_2(g1_tr_2_end+1:g1_tr_10_start-1,:); TS_COMPLETE_2(g2_tr_2_end+1:g2_tr_10_start-1,:); TS_COMPLETE_2(g3_tr_2_end+1:g3_tr_10_start-1,:); TS_COMPLETE_2(g4_tr_2_end+1:g4_tr_10_start-1,:); ];
L_TS_COMPLETE_2=[L_TS_COMPLETE_2(g1_tr_2_end+1:g1_tr_10_start-1,:); L_TS_COMPLETE_2(g2_tr_2_end+1:g2_tr_10_start-1,:); L_TS_COMPLETE_2(g3_tr_2_end+1:g3_tr_10_start-1,:); L_TS_COMPLETE_2(g4_tr_2_end+1:g4_tr_10_start-1,:); ];
% figure; hold on; 
% plot(TS_COMPLETE_2); plot(L_TS_COMPLETE_2)
% figure; hold on; 
% plot(SAMPL_DATA_2); plot(L_SAMPL_DATA_2)

%subject 3
g1_tr_2_end=33;
g1_tr_10_start=117;
g2_tr_2_end=166;
g2_tr_10_start=251;
g3_tr_2_end=314;
g3_tr_10_start=399;
g4_tr_2_end=446;
g4_tr_10_start=530;

SAMPL_DATA_3= [TS_COMPLETE_3(1:g1_tr_2_end,:);TS_COMPLETE_3(g1_tr_10_start:g2_tr_2_end,:); TS_COMPLETE_3(g2_tr_10_start:g3_tr_2_end,:); TS_COMPLETE_3(g3_tr_10_start:g4_tr_2_end,:); TS_COMPLETE_3(g4_tr_10_start:end,:);  ];
L_SAMPL_DATA_3= [L_TS_COMPLETE_3(1:g1_tr_2_end,:);L_TS_COMPLETE_3(g1_tr_10_start:g2_tr_2_end,:); L_TS_COMPLETE_3(g2_tr_10_start:g3_tr_2_end,:); L_TS_COMPLETE_3(g3_tr_10_start:g4_tr_2_end,:); L_TS_COMPLETE_3(g4_tr_10_start:end,:);];
%assign new values to TS_COMPLETE and L_TS_COMPLETE
TS_COMPLETE_3=[TS_COMPLETE_3(g1_tr_2_end+1:g1_tr_10_start-1,:); TS_COMPLETE_3(g2_tr_2_end+1:g2_tr_10_start-1,:); TS_COMPLETE_3(g3_tr_2_end+1:g3_tr_10_start-1,:); TS_COMPLETE_3(g4_tr_2_end+1:g4_tr_10_start-1,:); ];
L_TS_COMPLETE_3=[L_TS_COMPLETE_3(g1_tr_2_end+1:g1_tr_10_start-1,:); L_TS_COMPLETE_3(g2_tr_2_end+1:g2_tr_10_start-1,:); L_TS_COMPLETE_3(g3_tr_2_end+1:g3_tr_10_start-1,:); L_TS_COMPLETE_3(g4_tr_2_end+1:g4_tr_10_start-1,:); ];
% figure; hold on; 
% plot(TS_COMPLETE_3); plot(L_TS_COMPLETE_3)
%  figure; hold on; 
% plot(SAMPL_DATA_3); plot(L_SAMPL_DATA_3)


%subject 4
g1_tr_2_end=36;
g1_tr_10_start=120;
g2_tr_2_end=173;
g2_tr_10_start=258;
g3_tr_2_end=307;
g3_tr_10_start=392;
g4_tr_2_end=436;
g4_tr_10_start=521;

SAMPL_DATA_4= [TS_COMPLETE_4(1:g1_tr_2_end,:);TS_COMPLETE_4(g1_tr_10_start:g2_tr_2_end,:); TS_COMPLETE_4(g2_tr_10_start:g3_tr_2_end,:); TS_COMPLETE_4(g3_tr_10_start:g4_tr_2_end,:); TS_COMPLETE_4(g4_tr_10_start:end,:);  ];
L_SAMPL_DATA_4= [L_TS_COMPLETE_4(1:g1_tr_2_end,:);L_TS_COMPLETE_4(g1_tr_10_start:g2_tr_2_end,:); L_TS_COMPLETE_4(g2_tr_10_start:g3_tr_2_end,:); L_TS_COMPLETE_4(g3_tr_10_start:g4_tr_2_end,:); L_TS_COMPLETE_4(g4_tr_10_start:end,:);];
%assign new values to TS_COMPLETE and L_TS_COMPLETE
TS_COMPLETE_4=[TS_COMPLETE_4(g1_tr_2_end+1:g1_tr_10_start-1,:); TS_COMPLETE_4(g2_tr_2_end+1:g2_tr_10_start-1,:); TS_COMPLETE_4(g3_tr_2_end+1:g3_tr_10_start-1,:); TS_COMPLETE_4(g4_tr_2_end+1:g4_tr_10_start-1,:); ];
L_TS_COMPLETE_4=[L_TS_COMPLETE_4(g1_tr_2_end+1:g1_tr_10_start-1,:); L_TS_COMPLETE_4(g2_tr_2_end+1:g2_tr_10_start-1,:); L_TS_COMPLETE_4(g3_tr_2_end+1:g3_tr_10_start-1,:); L_TS_COMPLETE_4(g4_tr_2_end+1:g4_tr_10_start-1,:); ];
% figure; hold on; 
% plot(TS_COMPLETE_4); plot(L_TS_COMPLETE_4)
%  figure; hold on; 
% plot(SAMPL_DATA_4); plot(L_SAMPL_DATA_4)

%subject 5
g1_tr_2_end=177;
g1_tr_10_start=386;
g2_tr_2_end=531;
g2_tr_10_start=840;
g3_tr_2_end=981;
g3_tr_10_start=1335;
g4_tr_2_end=1494;
g4_tr_10_start=1865;

SAMPL_DATA_5= [TS_COMPLETE_5(1:g1_tr_2_end,:);TS_COMPLETE_5(g1_tr_10_start:g2_tr_2_end,:); TS_COMPLETE_5(g2_tr_10_start:g3_tr_2_end,:); TS_COMPLETE_5(g3_tr_10_start:g4_tr_2_end,:); TS_COMPLETE_5(g4_tr_10_start:end,:);  ];
L_SAMPL_DATA_5= [L_TS_COMPLETE_5(1:g1_tr_2_end,:);L_TS_COMPLETE_5(g1_tr_10_start:g2_tr_2_end,:); L_TS_COMPLETE_5(g2_tr_10_start:g3_tr_2_end,:); L_TS_COMPLETE_5(g3_tr_10_start:g4_tr_2_end,:); L_TS_COMPLETE_5(g4_tr_10_start:end,:);];
%assign new values to TS_COMPLETE and L_TS_COMPLETE
TS_COMPLETE_5=[TS_COMPLETE_5(g1_tr_2_end+1:g1_tr_10_start-1,:); TS_COMPLETE_5(g2_tr_2_end+1:g2_tr_10_start-1,:); TS_COMPLETE_5(g3_tr_2_end+1:g3_tr_10_start-1,:); TS_COMPLETE_5(g4_tr_2_end+1:g4_tr_10_start-1,:); ];
L_TS_COMPLETE_5=[L_TS_COMPLETE_5(g1_tr_2_end+1:g1_tr_10_start-1,:); L_TS_COMPLETE_5(g2_tr_2_end+1:g2_tr_10_start-1,:); L_TS_COMPLETE_5(g3_tr_2_end+1:g3_tr_10_start-1,:); L_TS_COMPLETE_5(g4_tr_2_end+1:g4_tr_10_start-1,:); ];
% figure; hold on; 
% plot(TS_COMPLETE_5); plot(L_TS_COMPLETE_5)
%  figure; hold on; 
% plot(SAMPL_DATA_5); plot(L_SAMPL_DATA_5)

% [L_SAMPL_DATA_1, SAMPL_DATA_1] = genTrainData (TS_COMPLETE_1, L_TS_COMPLETE_1, learningFrac, 'inorder');
% [L_SAMPL_DATA_2, SAMPL_DATA_2] = genTrainData (TS_COMPLETE_2, L_TS_COMPLETE_2, learningFrac, 'inorder');
% [L_SAMPL_DATA_3, SAMPL_DATA_3] = genTrainData (TS_COMPLETE_3, L_TS_COMPLETE_3, learningFrac, 'inorder');
% [L_SAMPL_DATA_4, SAMPL_DATA_4] = genTrainData (TS_COMPLETE_4, L_TS_COMPLETE_4, learningFrac, 'inorder');
% [L_SAMPL_DATA_5, SAMPL_DATA_5] = genTrainData (TS_COMPLETE_5, L_TS_COMPLETE_5, learningFrac, 'inorder');


    Simulation =50; %number of simulation runs
    DENS=[0.01:0.01:0.5];
    
    
    %DK one loop is to run several simulation for different seed, another
    %loop is for different densities for biplar and binary encoding
    ACC_bip=zeros(5,length(DENS),Simulation);
    ACC_bin=zeros(5,length(DENS),Simulation);
    ACC_spr=zeros(5,length(DENS),Simulation);
    ACC_bipd=zeros(5,length(DENS),Simulation);
    
 %%    
    for i=1:Simulation %
        seed=i;
        for j=1:length(DENS);
            density=DENS(j);
            disp([seed,density])

    codingMode = 'dense_binary';
    trainingMode = 'binary_vect';
    bipolarDistretizationMode = 'off';
    %iM = initItemMemories (D, 21, 6,seed,density); %DK. do iM once to save time for the largest N
    iM = initItemMemoriesLinear(D, 21, 6,seed,density);
    
    N = 4; %N-gram size 
    %Subject 1
    %iM = initItemMemories (D, 21, N);
    %iM = initItemMemories (D, 21, N,seed,density);
	[numpat, hdc_model_1] = hdctrain (L_SAMPL_DATA_1, SAMPL_DATA_1, iM, D, N, percision, cuttingAngle);
    [accExcTrnz, acch] = hdcpredict (L_TS_COMPLETE_1, TS_COMPLETE_1, hdc_model_1, iM, D, N, percision);
	%acc_holo_LF_1 (LF) = acch;
	%numpat_n_1 (N,:) = [numpat sum(numpat)];
	ACC_bin(1,j,i)=acch;
    
    
    N = 3; %N-gram size 
    %Subject 2
    %iM = initItemMemories (D, 21, N);
    %iM = initItemMemories (D, 21, N,seed,density);
    [numpat, hdc_model_2] = hdctrain (L_SAMPL_DATA_2, SAMPL_DATA_2, iM, D, N, percision, cuttingAngle);
    [accExcTrnz, acch] = hdcpredict (L_TS_COMPLETE_2, TS_COMPLETE_2, hdc_model_2, iM, D, N, percision);
	%acc_holo_LF_2 (LF) = acch;
    ACC_bin(2,j,i)=acch;
	
    %N = 7;
    N = 3; %N-gram size 
    %Subject 3
    %iM = initItemMemories (D, 21, N);
    %iM = initItemMemories (D, 21, N,seed,density);
    [numpat, hdc_model_3] = hdctrain (L_SAMPL_DATA_3, SAMPL_DATA_3, iM, D, N, percision, cuttingAngle);
    [accExcTrnz, acch] = hdcpredict (L_TS_COMPLETE_3, TS_COMPLETE_3, hdc_model_3, iM, D, N, percision);
	%acc_holo_LF_3 (LF) = acch;
    ACC_bin(3,j,i)=acch;
	
    N = 6; %N-gram size 
    %Subject 4
    %iM = initItemMemories (D, 21, N);
    %iM = initItemMemories (D, 21, N,seed,density);
    [numpat, hdc_model_4] = hdctrain (L_SAMPL_DATA_4, SAMPL_DATA_4, iM, D, N, percision, cuttingAngle);
    [accExcTrnz, acch] = hdcpredict (L_TS_COMPLETE_4, TS_COMPLETE_4, hdc_model_4, iM, D, N, percision);
	%acc_holo_LF_4 (LF) = acch;
    ACC_bin(4,j,i)=acch;
	
    N = 5; %N-gram size 
    %Subject 5
    %iM = initItemMemories (D, 21, N);
   % iM = initItemMemories (D, 21, N,seed,density);
    [numpat, hdc_model_5] = hdctrain (L_SAMPL_DATA_5, SAMPL_DATA_5, iM, D, N, percision, cuttingAngle);
    [accExcTrnz, acch] = hdcpredict (L_TS_COMPLETE_5, TS_COMPLETE_5, hdc_model_5, iM, D, N, percision);
	%acc_holo_LF_5 (LF) = acch;
    ACC_bin(5,j,i)=acch;
    %% 
    
    %DK. DO the same but for bipolar vectors
    codingMode = 'dense_bipolar';
    trainingMode = 'integer_vect';
    bipolarDistretizationMode = 'off';
        %iM = initItemMemories (D, 21, 6,seed,density); %DK. do iM once to save time for the largest N
        iM = initItemMemoriesLinear(D, 21, 6,seed,density);
        
        N = 4; %N-gram size 
    %Subject 1
    %iM = initItemMemories (D, 21, N);
    %iM = initItemMemories (D, 21, N,seed,density);
	[numpat, hdc_model_1] = hdctrain (L_SAMPL_DATA_1, SAMPL_DATA_1, iM, D, N, percision, cuttingAngle);
    [accExcTrnz, acch] = hdcpredict (L_TS_COMPLETE_1, TS_COMPLETE_1, hdc_model_1, iM, D, N, percision);
	%acc_holo_LF_1 (LF) = acch;
	%numpat_n_1 (N,:) = [numpat sum(numpat)];
	ACC_bip(1,j,i)=acch;
    
    
    N = 3; %N-gram size 
    %Subject 2
    %iM = initItemMemories (D, 21, N);
    %iM = initItemMemories (D, 21, N,seed,density);
    [numpat, hdc_model_2] = hdctrain (L_SAMPL_DATA_2, SAMPL_DATA_2, iM, D, N, percision, cuttingAngle);
    [accExcTrnz, acch] = hdcpredict (L_TS_COMPLETE_2, TS_COMPLETE_2, hdc_model_2, iM, D, N, percision);
	%acc_holo_LF_2 (LF) = acch;
    ACC_bip(2,j,i)=acch;
	
    %N = 7;
    N = 3; %N-gram size 
    %Subject 3
    %iM = initItemMemories (D, 21, N);
    %iM = initItemMemories (D, 21, N,seed,density);
    [numpat, hdc_model_3] = hdctrain (L_SAMPL_DATA_3, SAMPL_DATA_3, iM, D, N, percision, cuttingAngle);
    [accExcTrnz, acch] = hdcpredict (L_TS_COMPLETE_3, TS_COMPLETE_3, hdc_model_3, iM, D, N, percision);
	%acc_holo_LF_3 (LF) = acch;
    ACC_bip(3,j,i)=acch;
	
    N = 6; %N-gram size 
    %Subject 4
    %iM = initItemMemories (D, 21, N);
    %iM = initItemMemories (D, 21, N,seed,density);
    [numpat, hdc_model_4] = hdctrain (L_SAMPL_DATA_4, SAMPL_DATA_4, iM, D, N, percision, cuttingAngle);
    [accExcTrnz, acch] = hdcpredict (L_TS_COMPLETE_4, TS_COMPLETE_4, hdc_model_4, iM, D, N, percision);
	%acc_holo_LF_4 (LF) = acch;
    ACC_bip(4,j,i)=acch;
	
    N = 5; %N-gram size 
    %Subject 5
    %iM = initItemMemories (D, 21, N);
    %iM = initItemMemories (D, 21, N,seed,density);
    [numpat, hdc_model_5] = hdctrain (L_SAMPL_DATA_5, SAMPL_DATA_5, iM, D, N, percision, cuttingAngle);
    [accExcTrnz, acch] = hdcpredict (L_TS_COMPLETE_5, TS_COMPLETE_5, hdc_model_5, iM, D, N, percision);
	%acc_holo_LF_5 (LF) = acch;
    ACC_bip(5,j,i)=acch;

%% 
   
    %DK. DO the same but for bipolarized bipolar vectors
    codingMode = 'dense_bipolar';
    trainingMode = 'integer_vect';
    bipolarDistretizationMode = 'on';
    
    %global testngrams;
    
        %iM = initItemMemories (D, 21, 6,seed,density); %DK. do iM once to save time for the largest N
        iM = initItemMemoriesLinear(D, 21, 6,seed,density);
        
        N = 4; %N-gram size 
    %Subject 1
    %iM = initItemMemories (D, 21, N);
    %iM = initItemMemories (D, 21, N,seed,density);
	[numpat, hdc_model_1] = hdctrain (L_SAMPL_DATA_1, SAMPL_DATA_1, iM, D, N, percision, cuttingAngle);
    [accExcTrnz, acch] = hdcpredict (L_TS_COMPLETE_1, TS_COMPLETE_1, hdc_model_1, iM, D, N, percision);
	%acc_holo_LF_1 (LF) = acch;
	%numpat_n_1 (N,:) = [numpat sum(numpat)];
	ACC_bipd(1,j,i)=acch;
    
    
    N = 3; %N-gram size 
    %Subject 2
    %iM = initItemMemories (D, 21, N);
    %iM = initItemMemories (D, 21, N,seed,density);
    [numpat, hdc_model_2] = hdctrain (L_SAMPL_DATA_2, SAMPL_DATA_2, iM, D, N, percision, cuttingAngle);
    [accExcTrnz, acch] = hdcpredict (L_TS_COMPLETE_2, TS_COMPLETE_2, hdc_model_2, iM, D, N, percision);
	%acc_holo_LF_2 (LF) = acch;
    ACC_bipd(2,j,i)=acch;
	
    %N = 7;
    N = 3; %N-gram size 
    %Subject 3
    %iM = initItemMemories (D, 21, N);
    %iM = initItemMemories (D, 21, N,seed,density);
    [numpat, hdc_model_3] = hdctrain (L_SAMPL_DATA_3, SAMPL_DATA_3, iM, D, N, percision, cuttingAngle);
    [accExcTrnz, acch] = hdcpredict (L_TS_COMPLETE_3, TS_COMPLETE_3, hdc_model_3, iM, D, N, percision);
	%acc_holo_LF_3 (LF) = acch;
    ACC_bipd(3,j,i)=acch;
	
    N = 6; %N-gram size 
    %Subject 4
    %iM = initItemMemories (D, 21, N);
    %iM = initItemMemories (D, 21, N,seed,density);
    [numpat, hdc_model_4] = hdctrain (L_SAMPL_DATA_4, SAMPL_DATA_4, iM, D, N, percision, cuttingAngle);
    [accExcTrnz, acch] = hdcpredict (L_TS_COMPLETE_4, TS_COMPLETE_4, hdc_model_4, iM, D, N, percision);
	%acc_holo_LF_4 (LF) = acch;
    ACC_bipd(4,j,i)=acch;
	
    N = 5; %N-gram size 
    %Subject 5
    %iM = initItemMemories (D, 21, N);
    %iM = initItemMemories (D, 21, N,seed,density);
    [numpat, hdc_model_5] = hdctrain (L_SAMPL_DATA_5, SAMPL_DATA_5, iM, D, N, percision, cuttingAngle);
    [accExcTrnz, acch] = hdcpredict (L_TS_COMPLETE_5, TS_COMPLETE_5, hdc_model_5, iM, D, N, percision);
	%acc_holo_LF_5 (LF) = acch;
    ACC_bipd(5,j,i)=acch;


    %% Sparse coder
    N = 4; %N-gram size 
    %Subject 1
    
    codingMode = 'dense_binary';
    trainingMode = 'binary_vect';
    bipolarDistretizationMode = 'off';
    
    %iM = initItemMemoriesSparse (D, 21, 6, density, seed);
    iM = initItemMemoriesLinear(D, 21, 6,seed,density);
    
%iM = initItemMemoriesSparse (D, 21, N,100,i);
[numpat, hdc_model_1] = hdctrainSparse (L_SAMPL_DATA_1, SAMPL_DATA_1, iM, D, N, percision, cuttingAngle);
[accExcTrnz, acch_1] = hdcpredictSparse(L_TS_COMPLETE_1, TS_COMPLETE_1, hdc_model_1, iM, D, N, percision);

ACC_spr(1,j,i)=acch_1; 
%STAT_aet_spr(1,i)=accExcTrnz;
%STAT_acc_spr(1,i)=acch_1;

%[accExcTrnz, acch_1,probab1,PR] = hdcpredictconf(L_TS_COMPLETE_1, TS_COMPLETE_1, hdc_model_1, iM, D, N, percision);
%[accExcTrnz, acch_1,probab1,PR] = hdcpredictconf(LABEL_1, COMPLETE_1, hdc_model_1, iM, D, N, percision);


N = 3; %N-gram size 
    %Subject 2

%iM = initItemMemoriesSparse (D, 21, N,100,i);
[numpat, hdc_model_2] = hdctrainSparse (L_SAMPL_DATA_2, SAMPL_DATA_2, iM, D, N, percision, cuttingAngle);
[accExcTrnz, acch_2] = hdcpredictSparse(L_TS_COMPLETE_2, TS_COMPLETE_2, hdc_model_2, iM, D, N, percision);

ACC_spr(2,j,i)=acch_2;
%STAT_aet_spr(2,i)=accExcTrnz;
%STAT_acc_spr(2,i)=acch_2;

%[accExcTrnz, acch_2,probab2] = hdcpredictconf(L_TS_COMPLETE_2, TS_COMPLETE_2, hdc_model_2, iM, D, N, percision);

N = 3; %N-gram size 
    %Subject 3

%iM = initItemMemoriesSparse (D, 21, N,100,i);
[numpat, hdc_model_3] = hdctrainSparse (L_SAMPL_DATA_3, SAMPL_DATA_3, iM, D, N, percision, cuttingAngle);
[accExcTrnz, acch_3] = hdcpredictSparse(L_TS_COMPLETE_3, TS_COMPLETE_3, hdc_model_3, iM, D, N, percision);

ACC_spr(3,j,i)=acch_3;
%STAT_aet_spr(3,i)=accExcTrnz;
%STAT_acc_spr(3,i)=acch_3;


%[accExcTrnz, acch_3,probab3] = hdcpredictconf(L_TS_COMPLETE_3, TS_COMPLETE_3, hdc_model_3, iM, D, N, percision);

N = 6; %N-gram size 
    %Subject 4

%iM = initItemMemoriesSparse (D, 21, N,100,i);
[numpat, hdc_model_4] = hdctrainSparse (L_SAMPL_DATA_4, SAMPL_DATA_4, iM, D, N, percision, cuttingAngle);
[accExcTrnz, acch_4] = hdcpredictSparse(L_TS_COMPLETE_4, TS_COMPLETE_4, hdc_model_4, iM, D, N, percision);

ACC_spr(4,j,i)=acch_4;
%STAT_aet_spr(4,i)=accExcTrnz;
%STAT_acc_spr(4,i)=acch_4;

%[accExcTrnz, acch_4,probab4] = hdcpredictconf(L_TS_COMPLETE_4, TS_COMPLETE_4, hdc_model_4, iM, D, N, percision);

N = 5; %N-gram size 
    %Subject 5

%iM = initItemMemoriesSparse (D, 21, N,100,i);
[numpat, hdc_model_5] = hdctrainSparse (L_SAMPL_DATA_5, SAMPL_DATA_5, iM, D, N, percision, cuttingAngle);
[accExcTrnz, acch_5] = hdcpredictSparse(L_TS_COMPLETE_5, TS_COMPLETE_5, hdc_model_5, iM, D, N, percision);

ACC_spr(5,j,i)=acch_5;
%STAT_aet_spr(5,i)=accExcTrnz;
%STAT_acc_spr(5,i)=acch_5;  
 

        end
    
    end
    
%end
 
%for LF = st:res:sp    
%    acc_holo_LF_mean(LF) = mean([acc_holo_LF_1(LF) acc_holo_LF_2(LF) acc_holo_LF_3(LF) acc_holo_LF_4(LF) acc_holo_LF_5(LF)]);  
%    acc_holo_LF_std(LF) = std  ([acc_holo_LF_1(LF) acc_holo_LF_2(LF) acc_holo_LF_3(LF) acc_holo_LF_4(LF) acc_holo_LF_5(LF)]);  
%end


