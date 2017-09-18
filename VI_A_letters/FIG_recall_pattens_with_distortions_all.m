%This script outputs the figure corresponding to Fig. 5 in the article

warning off
clear all

%Provides set of images of letters. Fig. 2 in the original paper;
load Letters 

%Set the dimensionality of HD-vectors
d=10000;

d2=100000;
M2=1000;

d3=10000;
M3=100;

d4=2048;
M4=40;

%Set number of simulations for each noise level.
sim=1000;

%Number of GNs is determined as the number of pixels in the image
numGN=size(Letters{1,1}(:),1);

%Initialize arrays to deploy results of accuracy calculation
Accuracy_dense=zeros(5,26);
Accuracy_low=zeros(5,26);
Accuracy_spr_2k=zeros(5,26);
Accuracy_spr_10k=zeros(5,26);
Accuracy_spr_100k=zeros(5,26);
Noisy_all=cell(5,5);



%b is the number of bit errors. Range from 1 to maxb
maxb=5; %maximal number of bit errors

for b=1:maxb
   
   %Creates vector with b errors
   Err=[ones(b,1); zeros(numGN-b,1)];
          
   %Initialize array to deploy results of recalls
   Noisy_dense=zeros(sim,26); 
   Noisy_low=zeros(sim,26); 
   Noisy_spr_2k=zeros(sim,26);
   Noisy_spr_10k=zeros(sim,26);
   Noisy_spr_100k=zeros(sim,26);

   for k=1:sim
        disp([b,k]); 
       %Creates distributed representation (HGN) for every undistorted image of
        %letter via call to "letters_encoding" function and
        %"sparse_letters_encoding" function

        [IV_spr_2] = sparseVectorsGenerator(numGN,d2,M2);
        HGN_spr_2=sparse_letters_encoding(IV_spr_2,d2);
        [IV_spr_3] = sparseVectorsGenerator(numGN,d3,M3);
        HGN_spr_3=sparse_letters_encoding(IV_spr_3,d3);
        [IV_spr_4] = sparseVectorsGenerator(numGN,d4,M4);
        HGN_spr_4=sparse_letters_encoding(IV_spr_4,d4);       
        
        [IV_dn] = denseVectorsGenerator(numGN,d);
        [HGN, LD]=letters_encoding( d, IV_dn );
        %Complex array for later recall
        HGNc=bin2com(HGN);
        LDc=bin2com(LD);
        
        %Encodes every image of letter
        for j=1:26
        %Display current step in simulation
        %disp([b,j]); 
            rng('default');
            rng('shuffle')
            %Takes jth image from Letters  
            pict=Letters{j,1};

            %Reshapes image into pattern
            pattern(1,:)=pict(:)';   

            %Perform sim simulations for each noise level
       
            %Randomize error pattern
            Errn=Err(randperm(numGN));
            
            %Introduce distorions into image
            pattern_nois=double(xor(pattern,Errn'));
                        
            %Creates distributed representation for a noisy pattern
            repr=hologn_encoder(pattern_nois,d,IV_dn)';
            repr_spr_2=sparse_hologn_encoder_IV(pattern_nois,IV_spr_2)';
            repr_spr_3=sparse_hologn_encoder_IV(pattern_nois,IV_spr_3)';
            repr_spr_4=sparse_hologn_encoder_IV(pattern_nois,IV_spr_4)';
            
            
            %Recall the closest letter from the item memory which contains
            %distributed represetnations for undistorted letters (HGNc here)
            [Noisy_dense(k,j)]=item_memory_c(bin2com(repr), HGNc );
            [Noisy_low(k,j) ]=item_memory_c(bin2com(pattern_nois'), LDc ); 
            Noisy_spr_100k(k,j)=item_memory_overlap(repr_spr_2, HGN_spr_2 );
            Noisy_spr_10k(k,j) =item_memory_overlap(repr_spr_3, HGN_spr_3 );
            Noisy_spr_2k(k,j)  =item_memory_overlap(repr_spr_4, HGN_spr_4 );
            
        
       
        end    
   end
   Noisy_all{b,1}=Noisy_low;
   Noisy_all{b,2}=Noisy_dense;
   Noisy_all{b,3}=Noisy_spr_100k;
   Noisy_all{b,4}=Noisy_spr_10k;
   Noisy_all{b,5}=Noisy_spr_2k;
   
   %Calculate the percentage of average accuracy for noisy recall of every
   %letter for every level of distorsion
   %Note that every plot in Fig. 6. in the paper displays Accuracy(i,:), 
   %where i is the number of distorted bits  
    for i=1:26
        Accuracy_dense(b,i)=(sum((Noisy_dense(:,i)==i))/sim)*100;
        Accuracy_low(b,i)=(sum((Noisy_low(:,i)==i))/sim)*100; 
        Accuracy_spr_100k(b,i)=(sum((Noisy_spr_100k(:,i)==i))/sim)*100; 
        Accuracy_spr_10k(b,i) =(sum((Noisy_spr_10k(:,i)==i))/sim)*100; 
        Accuracy_spr_2k(b,i)  =(sum((Noisy_spr_2k(:,i)==i))/sim)*100; 
        
    end   
   
end



%% Plot histograms for each distortion value 

figure()
count=1;    
for b=1:2:maxb
    subplot(2,2,count)
    count=count+1;
    hold on
    
    coef=1;   
    x_ax=coef*1:coef*1:26*coef; % x axis
    bar( [x_ax]'   ,[Accuracy_low(b,:)/100; Accuracy_dense(b,:)/100; Accuracy_spr_100k(b,:)/100; Accuracy_spr_10k(b,:)/100; Accuracy_spr_2k(b,:)/100;   ]') 
    
    
    
    %plot confidence intervals    
    for j=1:5
        
        % dependent on the bar choose x-axis
        switch j
            case 1
                X=0.694*coef:coef*1:26.6*coef ;
            case 2
                X=0.846:26.6 ;
            case 3
                X=1.0:26.6 ;
            case 4
                X=1.153:26.6 ;
            case 5
                X=1.306:26.6 ;
        end        
        
        
        
        ACC=zeros(sim,26);
        for i=1:26
            ACC(:,i)=double(Noisy_all{b,j}(:,i)==i);        
        end
        ACC_st=std(ACC);
        ACC_mean=mean(ACC);
        interval= tinv(0.975,sim-1)*ACC_st./sqrt(sim);
        for i=1:26
            %errorbar(X,mean(ACC),interval,interval,'k.', 'LineWidth',2); 
            plot([X(i),X(i)],([ACC_mean(i)-interval(i), ACC_mean(i)+interval(i)]),'r','LineWidth',2 )          
        end
        
        %errorbar(X,mean(ACC),interval,interval,'k.', 'LineWidth',2);   
    end
    
    
    ax = gca;
    ax.XTick = (x_ax);
    ax.XTickLabel = ({'A','B','C','D', 'E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'});
    
    xlim([x_ax(1)-coef/2,x_ax(end)+coef/2]) 
    ylim([0,1]) 
    box on
    set(gca,'fontsize',20)
    xlabel('Letters','FontSize', 22)
    %ylabel('Result, %')
    ylabel('Accuracy','FontSize', 22)
    
    switch b
        case 1
            title('Distortion 1 bit (2.9%)')
        case 2
            title('Distortion 2 bits (5.7%)')
        case 3
            title('Distortion 3 bits (8.6%)')
        case 4
            title('Distortion 4 bits (11.4%)')
        case 5
            title('Distortion 5 bits (14.3%)')
            legend('low-dim representation', 'dense, 10 000 elements', 'sparse, 100 000 elements', 'sparse, 10 000 elements', 'sparse, 2048 elements')
    end
    
   
end






