%This script outputs the figure corresponding to Fig. 8 in the article

clear all
compareHoloGN %does the experiments for bipolar and binary with ideallinear mapping

n=50; %number of simulations
DENS=[0.01:0.01:0.5]; % density of ones used in simulations

%bipolar; linear; MAP; integer
ACC_bip_m=mean(ACC_bip,3); % estimation of the mean
ACC_bip_std=std(ACC_bip,0,3); % estimation of the standard deviation
ACC_bip_int= tinv(0.975,n-1)*ACC_bip_std/sqrt(n); %In order to estimate confidence intervals with 95% one need a coefficient from student-t distribution and with n-1 degrees of freedom where n is the number of simulations %interval= tinv(0.975,n-1)*S/sqrt(n), S is the estimation of the standard deviation
ACC_bip_up=ACC_bip_m+ACC_bip_int; % upper interval
ACC_bip_down=ACC_bip_m-ACC_bip_int; %lower interval


compareHoloGN_bipolarized %does the experiments for bipolarized
%bipolar; linear; majority rule; bipolar
ACC_bipd_m=mean(ACC_bipd,3); % estimation of the mean
ACC_bipd_std=std(ACC_bipd,0,3); % estimation of the standard deviation
ACC_bipd_int= tinv(0.975,n-1)*ACC_bipd_std/sqrt(n); %In order to estimate confidence intervals with 95% one need a coefficient from student-t distribution and with n-1 degrees of freedom where n is the number of simulations %interval= tinv(0.975,n-1)*S/sqrt(n), S is the estimation of the standard deviation
ACC_bipd_up=ACC_bipd_m+ACC_bipd_int; % upper interval
ACC_bipd_down=ACC_bipd_m-ACC_bipd_int; %lower interval



%binary; linear; majority rule; binary
ACC_bin_m=mean(ACC_bin,3); % estimation of the mean
ACC_bin_std=std(ACC_bin,0,3); % estimation of the standard deviation
ACC_bin_int= tinv(0.975,n-1)*ACC_bin_std/sqrt(n); %In order to estimate confidence intervals with 95% one need a coefficient from student-t distribution and with n-1 degrees of freedom where n is the number of simulations %interval= tinv(0.975,n-1)*S/sqrt(n), S is the estimation of the standard deviation
ACC_bin_up=ACC_bin_m+ACC_bin_int; % upper interval
ACC_bin_down=ACC_bin_m-ACC_bin_int; %lower interval

compareHoloGN_sparse %does the experiments for sparse
%binary; linear; CDT; binary
ACC_spr_m=mean(ACC_spr,3); % estimation of the mean
ACC_spr_std=std(ACC_spr,0,3); % estimation of the standard deviation
ACC_spr_int= tinv(0.975,n-1)*ACC_spr_std/sqrt(n); %In order to estimate confidence intervals with 95% one need a coefficient from student-t distribution and with n-1 degrees of freedom where n is the number of simulations %interval= tinv(0.975,n-1)*S/sqrt(n), S is the estimation of the standard deviation
ACC_spr_up=ACC_spr_m+ACC_spr_int; % upper interval
ACC_spr_down=ACC_spr_m-ACC_spr_int; %lower interval




%Plot all results
figure()
for i=1:5

    
    subplot(2,3,i)
    hold on
    
    %plot confidence intervals
    
    hpatch=patch([DENS fliplr(DENS)], [ ACC_bip_down(i,:) fliplr(ACC_bip_up(i,:))], 0.50*[1 1 1]); %bipolar; linear; MAP; integer
    set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.5);    

    hpatch=patch([DENS fliplr(DENS)], [ ACC_bipd_down(i,:) fliplr(ACC_bipd_up(i,:))], 0.50*[1 1 1]); %bipolar; linear; majority rule; bipolar
    set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.5);       
    
    hpatch=patch([DENS fliplr(DENS)], [ ACC_bin_down(i,:) fliplr(ACC_bin_up(i,:))], 0.50*[1 1 1]); %binary; linear; majority rule; binary 
    set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    hpatch=patch([DENS fliplr(DENS)], [ ACC_spr_down(i,:) fliplr(ACC_spr_up(i,:))], 0.50*[1 1 1]);  %binary; linear; CDT; binary
    set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.5);    
   
    
    %plot mean values
    l1=plot(DENS, ACC_bip_m(i,:),'-.g','Linewidth',3); %bipolar; linear; MAP; integer
    l2=plot(DENS, ACC_bipd_m(i,:),'r','Linewidth',3); %bipolar; linear; majority rule; bipolar
    l3=plot(DENS, ACC_bin_m(i,:),'--b','Linewidth',3); %binary; linear; majority rule; binary 
    l4=plot(DENS, ACC_spr_m(i,:),':k','Linewidth',3); %binary; linear; CDT; binary
    
    

    set(gca,'fontsize',20)
    box on
    grid on
    xlim([0,0.5])
    ylim([0.4,1.0])
    xlabel('Density of ones','FontSize', 22)
    ylabel('Accuracy','FontSize', 22)
    
    switch i
            case 1
                title('Subject 1','FontSize', 22)
            case 2
                title('Subject 2')
            case 3
                title('Subject 3')
            case 4
                title('Subject 4')
            case 5
                title('Subject 5')
                legend([l1 l2 l3 l4],{'bipolar; linear; MAP; integer','bipolar; linear; majority rule; bipolar','binary; linear; majority rule; binary ','binary; linear; CDT; binary'}); 
                %legend('show')
                %legend('low-dim representation', 'dense, 10 000 elements', 'sparse, 100 000 elements', 'sparse, 10 000 elements', 'sparse, 2048 elements')
    end
    
    
    
    
    
end

