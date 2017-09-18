%This script outputs the figure corresponding to Fig. 7 in the article

    clear all

    %Provides set of images of letters. Fig. 4 in the original paper;
    load Letters 

    %Set the dimensionality of HD-vectors
    d=10000;
    
    d2=100000;
    M2=1000; %Expected number of nonzero elements in an HD vector

    d3=10000;
    M3=100; %Expected number of nonzero elements in an HD vector

    d4=2048;
    M4=40; %Expected number of nonzero elements in an HD vector
    
    %Number of bit errors
    b=5; 
    
    %Creates vector with b errors
    Err=[ones(b,1); zeros(35-b,1)];
    
    %Number of GNs is determined as the number of pixels in the image
    numGN=size(Letters{1,1}(:),1);
    
    %Set the number of recalls
    Noisy_recall=100;
    
    %Set the maximal size of the training  set
    Bundle_limit=50;
    
    % Initialize cell array for storage of the training set
    BUNDLE=cell(26,5);
    
    %Initialize array to deploy results of simulation
    Accuracy_Noisy_ds=zeros(3,Bundle_limit);
    Accuracy_Noisy_ld=zeros(3,Bundle_limit);
    Accuracy_Noisy_spr_2k=zeros(3,Bundle_limit);
    Accuracy_Noisy_spr_10k=zeros(3,Bundle_limit);
    Accuracy_Noisy_spr_100k=zeros(3,Bundle_limit);
    
    
    %generate a dictionary of IV HD vectors
    [IV_dn] = denseVectorsGenerator(numGN,d);
    [IV_spr_2] = sparseVectorsGenerator(numGN,d2,M2);       
    [IV_spr_3] = sparseVectorsGenerator(numGN,d3,M3);       
    [IV_spr_4] = sparseVectorsGenerator(numGN,d4,M4);  
    Noisy_all=cell(Bundle_limit,5); % store all statistics
    
    %Do simulation for each size of the training set 
    for i=1:Bundle_limit
        %Display current step in simulation
        disp(i)
        % Initialize array for superposition of the training set
        HGN=zeros(26,d);
        LD=zeros(26,numGN);
        HGN_spr_2=zeros(26,d2);
        HGN_spr_3=zeros(26,d3);
        HGN_spr_4=zeros(26,d4);        
        
        
        %Initialize array to deploy results of recalls of noisy patterns 
        Noise_dense=zeros(Noisy_recall,26);
        Noise_ld=zeros(Noisy_recall,26);
        Noise_spr_2k=zeros(Noisy_recall,26);
        Noise_spr_10k=zeros(Noisy_recall,26);
        Noise_spr_100k=zeros(Noisy_recall,26);      
        
        %Encodes every image of letter
        for j=1:26 % number of letters
              rng('default');
              rng('shuffle');

              %Takes jth image from Letters  
              pict=Letters{j,1};

              %Reshapes image into pattern
              pattern=pict(:)';  

              %Randomize error pattern
              Errn=Err(randperm(numGN));

              %Introduce distorions into image
              pattern_nois=double(xor(pattern,Errn'));                

              %Creates distributed representation for a noisy pattern
              repr=hologn_encoder(pattern_nois,d,IV_dn)';              
              repr_spr_2=sparse_hologn_encoder_IV(pattern_nois,IV_spr_2)';
              repr_spr_3=sparse_hologn_encoder_IV(pattern_nois,IV_spr_3)';
              repr_spr_4=sparse_hologn_encoder_IV(pattern_nois,IV_spr_4)';
              
              %Add representation of a noisy pattetn into training set
              BUNDLE{j,1}(end+1,:)=repr; 
              BUNDLE{j,2}(end+1,:)=pattern_nois;
              BUNDLE{j,3}(end+1,:)=repr_spr_2;
              BUNDLE{j,4}(end+1,:)=repr_spr_3;
              BUNDLE{j,5}(end+1,:)=repr_spr_4;
              
              %Create superposition of training set through majority sum
              %operation
              HGN(j,:)=majority_sum(BUNDLE{j,1});
              LD(j,:)=(sum(BUNDLE{j,2})/i);
              HGN_spr_2(j,:)=Majority_CDT(BUNDLE{j,3},M2,d2);
              HGN_spr_3(j,:)=Majority_CDT(BUNDLE{j,4},M3,d3);
              HGN_spr_4(j,:)=Majority_CDT(BUNDLE{j,5},M4,d4);
        end

        %Perform recall of noisy pattern on trained representations 
        %recall is done Noisy_recall times
        parfor k=1:Noisy_recall
            %Display current step in simulation
            %disp([i,j,k]);
            
            %Encodes every image of letter
            for j=1:26 % number of letters
                rng('default');
                rng('shuffle');
                
                %Takes jth image from Letters                 
                pict=Letters{j,1};

                %Reshapes image into pattern
                pattern=pict(:)';  
                
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
                %distributed represetnations extracted from the training set (HGN here)
                Noise_dense(k,j)=item_memory_c(bin2com(repr), bin2com(HGN) );
                Noise_ld(k,j)=item_memory_ed(pattern_nois', LD );  
                Noise_spr_100k(k,j)=item_memory_overlap(repr_spr_2, HGN_spr_2 );
                Noise_spr_10k(k,j) =item_memory_overlap(repr_spr_3, HGN_spr_3 );
                Noise_spr_2k(k,j)  =item_memory_overlap(repr_spr_4, HGN_spr_4 );
                
            end            

        end
        
        %collect statistics
        Noisy_all(i,:)={Noise_ld, Noise_dense,Noise_spr_100k,Noise_spr_10k, Noise_spr_2k  };
        
        %Calculate the percentage of average accuracy for noisy recall of every
        %letter for every size of training set
        Succ_rate_ds=zeros(1,26);
        Succ_rate_ld=zeros(1,26);
        Succ_rate_spr_100k=zeros(1,26);
        Succ_rate_spr_10k=zeros(1,26);
        Succ_rate_spr_2k=zeros(1,26);
        
        for j=1:26
            Succ_rate_ds(1,j)=(sum((Noise_dense(:,j)==j))/Noisy_recall)*100;
            Succ_rate_ld(1,j)=(sum((Noise_ld(:,j)==j))/Noisy_recall)*100;
            Succ_rate_spr_100k(1,j)=(sum((Noise_spr_100k(:,j)==j))/Noisy_recall)*100;
            Succ_rate_spr_10k(1,j)=(sum((Noise_spr_10k(:,j)==j))/Noisy_recall)*100;
            Succ_rate_spr_2k(1,j)=(sum((Noise_spr_2k(:,j)==j))/Noisy_recall)*100;
        end
        
        %Calculate the percentage of average accuracy for noisy recall for
        %every size of training set along with maximum and minumum
        %accuracires between letters
        %Note that plot in Fig. 8. in the paper displays Accuracy_Noisy(1,:) 
        Accuracy_Noisy_ds(:,i)=[mean(Succ_rate_ds);max(Succ_rate_ds);min(Succ_rate_ds)] ;
        Accuracy_Noisy_ld(:,i)=[mean(Succ_rate_ld);max(Succ_rate_ld);min(Succ_rate_ld)] ;
        Accuracy_Noisy_spr_2k(:,i)=[mean(Succ_rate_spr_2k);max(Succ_rate_spr_2k);min(Succ_rate_spr_2k)] ;
        Accuracy_Noisy_spr_10k(:,i)=[mean(Succ_rate_spr_10k);max(Succ_rate_spr_10k);min(Succ_rate_spr_10k)] ;
        Accuracy_Noisy_spr_100k(:,i)=[mean(Succ_rate_spr_100k);max(Succ_rate_spr_100k);min(Succ_rate_spr_100k)] ;
        
    end
    
    

    
    

    
    %%     
    figure()
    hold on

    
    
     %Process collected statistics in order to get confidence intervals
    for j=1:5
        for i=1:Bundle_limit
            ACC=zeros(Noisy_recall,26);
            for k=1:26
                ACC(:,k)=double(Noisy_all{i,j}(:,k)==k);        
            end
            ACC=ACC(:);
            MEAN(i,j)=mean(ACC);
            CONF(i,j)= tinv(0.975,26*Noisy_recall-1)*std(ACC)./sqrt(26*Noisy_recall);   
        end
        
        
    end
    
    % depict confidence intervals as shaded areas
    DIST=1:Bundle_limit;
    j=2; % dense
    hpatch=patch([DIST fliplr(DIST)], [ (MEAN(:,j)-CONF(:,j))' fliplr((MEAN(:,j)+CONF(:,j))' )], 0.50*[0 0 0]); %bipolar; linear; MAP; integer
    set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.1);  
    j=5; % sparse 2,048
    hpatch=patch([DIST fliplr(DIST)], [ (MEAN(:,j)-CONF(:,j))' fliplr((MEAN(:,j)+CONF(:,j))' )], 0.50*[0 0 1]); %bipolar; linear; MAP; integer
    set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.1);  
    j=4; % sparse 10,000
    hpatch=patch([DIST fliplr(DIST)], [ (MEAN(:,j)-CONF(:,j))' fliplr((MEAN(:,j)+CONF(:,j))' )], 0.50*[1 0 0]); %bipolar; linear; MAP; integer
    set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.1);          
    j=3; % sparse 100,000
    hpatch=patch([DIST fliplr(DIST)], [ (MEAN(:,j)-CONF(:,j))' fliplr((MEAN(:,j)+CONF(:,j))' )], 0.50*[0 1 0]); %bipolar; linear; MAP; integer
    set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.1);      
    %plot confidence intervals 
    j=1; % low dimensional
    hpatch=patch([DIST fliplr(DIST)], [ (MEAN(:,j)-CONF(:,j))' fliplr((MEAN(:,j)+CONF(:,j))' )], 0.50*[1 0 1]); %bipolar; linear; MAP; integer
    set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.1);    
  
    l1=plot(Accuracy_Noisy_ld(1,:)/100,'m','LineWidth',2);  
    l2=plot(Accuracy_Noisy_ds(1,:)/100,'k','LineWidth',2);
    l3=plot(Accuracy_Noisy_spr_100k(1,:)/100,':g','LineWidth',2);
    l4=plot(Accuracy_Noisy_spr_10k(1,:)/100,'-.r','LineWidth',2);
    l5=plot(Accuracy_Noisy_spr_2k(1,:)/100,'--b','LineWidth',2);
    
  
    set(gca,'fontsize',20)
    xlim([0,50])
    ylim([0.4,1])
    ax = gca;
    ax.XTick = (0:5:50);
    
    xlabel('Size of training set','FontSize', 22)
    ylabel('Accuracy','FontSize', 22)
    grid on
    box on
    legend([l1 l2 l3 l4 l5],{'low-dim representation', 'dense, 10 000 elements', 'sparse, 100 000 elements', 'sparse, 10 000 elements', 'sparse, 2048 elements'});
    title('Distortion 5 bits (14.3%)')
        
    
 
    
    
    
    
    