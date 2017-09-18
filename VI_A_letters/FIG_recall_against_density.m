%This script outputs the figure corresponding to Fig. 4 in the article

figure
hold on

%% Calculate curves for dimensionality 100

warning off
clear all
d=100; %Set the dimensionality of HD-vectors
DENS=0.01:0.01:0.5;  % a range of density of ones
maxb=5; %range of bit distortions
sim=100; %Set number of simulations for each noise level.
Scenario_recall_pattens_with_distortions_all_density

%initialize arrays for storing statistics
Accuracy_dense=zeros(maxb,length(DENS)); 
Accuracy_sparse=zeros(maxb,length(DENS));
Interval_dense=zeros(maxb,length(DENS));
Interval_sparse=zeros(maxb,length(DENS));    

%% Plot curves for the chosen dimensionality 
for b=1:maxb
    

    
    for j=1:length(DENS)

        ACC_s=zeros(sim,26);
        ACC_d=zeros(sim,26);
        for i=1:26
             ACC_d(:,i)=double(Noisy_all_dense{b,j}(:,i)==i); 
             ACC_s(:,i)=double(Noisy_all_sparse{b,j}(:,i)==i);           
        end 
        ACC_d=mean(ACC_d,2);
        ACC_s=mean(ACC_s,2);
        ACC_d_st=std(ACC_d);
        ACC_s_st=std(ACC_s);
        %confidence interval
        Interval_dense(b,j)=tinv(0.975,sim-1)*ACC_d_st./sqrt(sim);
        Interval_sparse(b,j)=tinv(0.975,sim-1)*ACC_s_st./sqrt(sim);
        %mean accuracy
        Accuracy_dense(b,j)=mean(ACC_d);
        Accuracy_sparse(b,j)=mean(ACC_s);
    end
     

    %plot only for b=1;3;5
    if mod(b,2)==1
        
        
        
        subplot(1,3,(b+1)/2)
        hold on  
  
        %hpatch=patch([DENS fliplr(DENS)], [ (Accuracy_dense(b,:)-Interval_dense(b,:)) fliplr((Accuracy_dense(b,:)+Interval_dense(b,:)) )], 1.0*[1 1 1]); %bipolar; linear; MAP; integer
        %set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.1); 

        %hpatch=patch([DENS fliplr(DENS)], [ (Accuracy_sparse(b,:)-Interval_sparse(b,:)) fliplr((Accuracy_sparse(b,:)+Interval_sparse(b,:)) )], 1.0*[1 1 1]); %bipolar; linear; MAP; integer
        %set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.1); 

        plot(DENS,Accuracy_dense(b,:),'r','LineWidth',1.5)
        plot(DENS,Accuracy_sparse(b,:),'b','LineWidth',1.5)
        
        set(gca,'fontsize',20)
        xlim([0,0.5])
        ylim([0,1])
        xlabel('Density of ones','FontSize', 22)
        ylabel('Accuracy','FontSize', 22)
        box on    
        grid on
        
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
        end
        
        
    end
        
end
%% END Calculate curves for dimensionality 100

%% Calculate curves for dimensionality 1000

warning off
clear all
d=1000; %Set the dimensionality of HD-vectors
DENS=0.01:0.01:0.5;  % a range of density of ones
maxb=5; %range of bit distortions
sim=100; %Set number of simulations for each noise level.
Scenario_recall_pattens_with_distortions_all_density

%initialize arrays for storing statistics
Accuracy_dense=zeros(maxb,length(DENS)); 
Accuracy_sparse=zeros(maxb,length(DENS));
Interval_dense=zeros(maxb,length(DENS));
Interval_sparse=zeros(maxb,length(DENS));    

%% Plot curves for the chosen dimensionality 
for b=1:maxb
    

    
    for j=1:length(DENS)

        ACC_s=zeros(sim,26);
        ACC_d=zeros(sim,26);
        for i=1:26
             ACC_d(:,i)=double(Noisy_all_dense{b,j}(:,i)==i); 
             ACC_s(:,i)=double(Noisy_all_sparse{b,j}(:,i)==i);           
        end 
        ACC_d=mean(ACC_d,2);
        ACC_s=mean(ACC_s,2);
        ACC_d_st=std(ACC_d);
        ACC_s_st=std(ACC_s);
        %confidence interval
        Interval_dense(b,j)=tinv(0.975,sim-1)*ACC_d_st./sqrt(sim);
        Interval_sparse(b,j)=tinv(0.975,sim-1)*ACC_s_st./sqrt(sim);
        %mean accuracy
        Accuracy_dense(b,j)=mean(ACC_d);
        Accuracy_sparse(b,j)=mean(ACC_s);
    end
     

    %plot only for b=1;3;5
    if mod(b,2)==1
        
        
        
        subplot(1,3,(b+1)/2)
        hold on  
  
        %hpatch=patch([DENS fliplr(DENS)], [ (Accuracy_dense(b,:)-Interval_dense(b,:)) fliplr((Accuracy_dense(b,:)+Interval_dense(b,:)) )], 1.0*[1 1 1]); %bipolar; linear; MAP; integer
        %set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.1); 

        %hpatch=patch([DENS fliplr(DENS)], [ (Accuracy_sparse(b,:)-Interval_sparse(b,:)) fliplr((Accuracy_sparse(b,:)+Interval_sparse(b,:)) )], 1.0*[1 1 1]); %bipolar; linear; MAP; integer
        %set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.1); 

        plot(DENS,Accuracy_dense(b,:),'--r','LineWidth',1.5)
        plot(DENS,Accuracy_sparse(b,:),'--b','LineWidth',1.5)
        
        set(gca,'fontsize',20)
        xlim([0,0.5])
        ylim([0,1])
        xlabel('Density of ones','FontSize', 22)
        ylabel('Accuracy','FontSize', 22)
        box on    
        grid on
        
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
        end
        
        
    end
        
end
%% END Calculate curves for dimensionality 1000

%% Calculate curves for dimensionality 10, 000

warning off
clear all
d=10000; %Set the dimensionality of HD-vectors
DENS=0.01:0.01:0.5;  % a range of density of ones
maxb=5; %range of bit distortions
sim=100; %Set number of simulations for each noise level.
Scenario_recall_pattens_with_distortions_all_density

%initialize arrays for storing statistics
Accuracy_dense=zeros(maxb,length(DENS)); 
Accuracy_sparse=zeros(maxb,length(DENS));
Interval_dense=zeros(maxb,length(DENS));
Interval_sparse=zeros(maxb,length(DENS));    

%% Plot curves for the chosen dimensionality 
for b=1:maxb
    

    
    for j=1:length(DENS)

        ACC_s=zeros(sim,26);
        ACC_d=zeros(sim,26);
        for i=1:26
             ACC_d(:,i)=double(Noisy_all_dense{b,j}(:,i)==i); 
             ACC_s(:,i)=double(Noisy_all_sparse{b,j}(:,i)==i);           
        end 
        ACC_d=mean(ACC_d,2);
        ACC_s=mean(ACC_s,2);
        ACC_d_st=std(ACC_d);
        ACC_s_st=std(ACC_s);
        %confidence interval
        Interval_dense(b,j)=tinv(0.975,sim-1)*ACC_d_st./sqrt(sim);
        Interval_sparse(b,j)=tinv(0.975,sim-1)*ACC_s_st./sqrt(sim);
        %mean accuracy
        Accuracy_dense(b,j)=mean(ACC_d);
        Accuracy_sparse(b,j)=mean(ACC_s);
    end
     

    %plot only for b=1;3;5
    if mod(b,2)==1
        
        
        
        subplot(1,3,(b+1)/2)
        hold on  
  
        %hpatch=patch([DENS fliplr(DENS)], [ (Accuracy_dense(b,:)-Interval_dense(b,:)) fliplr((Accuracy_dense(b,:)+Interval_dense(b,:)) )], 1.0*[1 1 1]); %bipolar; linear; MAP; integer
        %set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.1); 

        %hpatch=patch([DENS fliplr(DENS)], [ (Accuracy_sparse(b,:)-Interval_sparse(b,:)) fliplr((Accuracy_sparse(b,:)+Interval_sparse(b,:)) )], 1.0*[1 1 1]); %bipolar; linear; MAP; integer
        %set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.1); 

        plot(DENS,Accuracy_dense(b,:),'-.r','LineWidth',1.5)
        plot(DENS,Accuracy_sparse(b,:),'-.b','LineWidth',1.5)
        
        set(gca,'fontsize',20)
        xlim([0,0.5])
        ylim([0,1])
        xlabel('Density of ones','FontSize', 22)
        ylabel('Accuracy','FontSize', 22)
        box on    
        grid on
        
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
        end
        
        
    end
        
end
%% END Calculate curves for dimensionality 10, 000