

%Provides set of images of letters. Fig. 2 in the original paper;
load Letters 

%Set the dimensionality of HD-vectors
%d=100000;

DENS=0.01:0.01:0.5; % a range of density of ones


%Set number of simulations for each noise level.
%sim=1000;

%Number of GNs is determined as the number of pixels in the image
numGN=size(Letters{1,1}(:),1);

%Initialize arrays to deploy results of accuracy calculation
% Accuracy_dense=zeros(5,26);
% Accuracy_low=zeros(5,26);
% Accuracy_spr_2k=zeros(5,26);
% Accuracy_spr_10k=zeros(5,26);
% Accuracy_spr_100k=zeros(5,26);

%b is the number of bit errors. Range from 1 to maxb
maxb=5; %maximal number of bit errors

Noisy_all_dense=cell(maxb,length(DENS)); %stores recall statistics
Noisy_all_sparse=cell(maxb,length(DENS)); %stores recall statistics

for j=1:length(DENS)
    rho=DENS(j);
    for b=1:maxb
       disp([j,b]);
       %Creates vector with b errors
       Err=[ones(b,1); zeros(numGN-b,1)];

       %Initialize array to deploy results of recalls
       Noisy_dense=zeros(sim,26);    
       Noisy_sparse=zeros(sim,26);


       for k=1:sim
            %disp([b,k]); 
           %Creates distributed representation (HGN) for every undistorted image of
            %letter via call to "letters_encoding" function and
            %"sparse_letters_encoding" function

            [IV] = sparseVectorsGenerator_sparsity(numGN,d,rho);
            
            HGN_spr=sparse_letters_encoding(IV,d);      
            [HGN]=letters_encoding( d, IV );
            %Complex array for later recall
            HGNc=bin2com(HGN);

            
            %Encodes every image of letter
            for i=1:26
            %Display current step in simulation
            %disp([b,j]); 
                rng('default');
                rng('shuffle')
                %Takes jth image from Letters  
                pict=Letters{i,1};

                %Reshapes image into pattern
                pattern=pict(:)';   

                %Perform sim simulations for each noise level

                %Randomize error pattern
                Errn=Err(randperm(numGN));

                %Introduce distorions into image
                pattern_nois=double(xor(pattern,Errn'));

                %Creates distributed representation for a noisy pattern
                repr_dense=hologn_encoder(pattern_nois,d,IV)';
                repr_sparse=sparse_hologn_encoder_IV(pattern_nois,IV)';


                %Recall the closest letter from the item memory which contains
                %distributed represetnations for undistorted letters (HGNc here)
                [Noisy_dense(k,i)]=item_memory_c(bin2com(repr_dense), HGNc ); 
                Noisy_sparse(k,i)=item_memory_overlap(repr_sparse, HGN_spr );

            end    
       end
       
       Noisy_all_dense{b,j}=Noisy_dense; %stores recall statistics
       Noisy_all_sparse{b,j}=Noisy_sparse; %stores recall statistics 
       
       
       %Calculate the percentage of average accuracy for noisy recall of every
       %letter for every level of distorsion
       %Note that every plot in Fig. 6. in the paper displays Accuracy(i,:), 
       %where i is the number of distorted bits  
%         for i=1:26
%             Accuracy_dense(b,i)=(sum((Noisy_dense(:,i)==i))/sim)*100;
%             Accuracy_low(b,i)=(sum((Noisy_low(:,i)==i))/sim)*100; 
%             Accuracy_spr_100k(b,i)=(sum((Noisy_spr_100k(:,i)==i))/sim)*100; 
%             Accuracy_spr_10k(b,i) =(sum((Noisy_spr_10k(:,i)==i))/sim)*100; 
%             Accuracy_spr_2k(b,i)  =(sum((Noisy_spr_2k(:,i)==i))/sim)*100; 
% 
%         end   

    end

end

%% Plot histograms for each distortion value 

% figure()
%     
% for b=1:maxb
%     subplot(3,2,b)
%     hold on
%     bar( [1:26]'   ,[Accuracy_low(b,:); Accuracy_dense(b,:); Accuracy_spr_100k(b,:); Accuracy_spr_10k(b,:); Accuracy_spr_2k(b,:);   ]') 
%     ax = gca;
%     ax.XTick = ([1:26]);
%     ax.XTickLabel = ({'A','B','C','D', 'E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'});
%     
%     xlim([0.5,26.5]) 
%     box on
%     xlabel('Letters')
%     ylabel('Result, %')
%     
%     switch b
%         case 1
%             title('Distortion 1 bit (2.9%)')
%         case 2
%             title('Distortion 2 bits (5.7%)')
%         case 3
%             title('Distortion 3 bits (8.6%)')
%         case 4
%             title('Distortion 4 bits (11.4%)')
%         case 5
%             title('Distortion 5 bits (14.3%)')
%             legend('low-dim representation', 'dense, 10 000 elements', 'sparse, 100 000 elements', 'sparse, 10 000 elements', 'sparse, 2048 elements')
%     end
%     
%    
% end




%In order to estimate confidence intervals with 95% one need a coefficient
%from student-t distribution and with n-1 degrees of freedom where n is the
%number of simulations


%interval= tinv(0.975,n-1)*S/sqrt(n), S is the estimation of the standard
%deviation




%% Plot histograms for each distortion value 

% figure(1)
%     
% for b=1:maxb
%     
%     Accuracy_dense=zeros(length(DENS),26);
%     Accuracy_sparse=zeros(length(DENS),26);
%     
%     for j=1:length(DENS)
% 
% 
%         for i=1:26
%              Accuracy_dense(j,i)=(sum((Noisy_all_dense{b,j}(:,i)==i))/sim); 
%              Accuracy_sparse(j,i)=(sum((Noisy_all_sparse{b,j}(:,i)==i))/sim);              
%         end 
%         
%         
%         
%         
%         %bar( [1:26]', [Accuracy_low(b,:); Accuracy_dense(b,:); Accuracy_spr_100k(b,:); Accuracy_spr_10k(b,:); Accuracy_spr_2k(b,:);   ]') 
%         %ax = gca;
%         %ax.XTick = ([1:26]);
%         %ax.XTickLabel = ({'A','B','C','D', 'E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'});
% 
%         %xlim([0.5,26.5]) 
% 
%     end
%     
%     
%     
%     subplot(3,2,b)
%     hold on
%     box on
%     plot(DENS,mean(Accuracy_dense,2),'-.r')
%     plot(DENS,mean(Accuracy_sparse,2),'-.b')
%     
%     xlabel('Sparsity')
%     ylabel('Result, %')
%          
%     
%     
%         switch b
%             case 1
%                 title('Distortion 1 bit (2.9%)')
%             case 2
%                 title('Distortion 2 bits (5.7%)')
%             case 3
%                 title('Distortion 3 bits (8.6%)')
%             case 4
%                 title('Distortion 4 bits (11.4%)')
%             case 5
%                 title('Distortion 5 bits (14.3%)')
%                 %legend('low-dim representation', 'dense, 10 000 elements', 'sparse, 100 000 elements', 'sparse, 10 000 elements', 'sparse, 2048 elements')
%         end
%         
% end

save SIMUL_NEW_d_100000_04_05
