%This script outputs the figure corresponding to Fig. 9 in the article


n=50; %number of simulations
DENS=[0.01:0.01:0.5]; % density of ones used in simulations



compareHoloGN %does the experiments for bipolar and binary with ideallinear mapping
%bipolar; linear; MAP; integer
ACC_bip_m=mean(ACC_bip,3); % estimation of the mean
ACC_bip_std=std(ACC_bip,0,3); % estimation of the standard deviation
ACC_bip_int= tinv(0.975,n-1)*ACC_bip_std/sqrt(n); %In order to estimate confidence intervals with 95% one need a coefficient from student-t distribution and with n-1 degrees of freedom where n is the number of simulations %interval= tinv(0.975,n-1)*S/sqrt(n), S is the estimation of the standard deviation
ACC_bip_up=ACC_bip_m+ACC_bip_int; % upper interval
ACC_bip_down=ACC_bip_m-ACC_bip_int; %lower interval



compareHoloGN_linear %does the experiments for approx. linear

%binary; approx. linear; majority rule; binary
ACC_bin_l_m=mean(ACC_bin,3); % estimation of the mean
ACC_bin_l_std=std(ACC_bin,0,3); % estimation of the standard deviation
ACC_bin_l_int= tinv(0.975,n-1)*ACC_bin_l_std/sqrt(n); %In order to estimate confidence intervals with 95% one need a coefficient from student-t distribution and with n-1 degrees of freedom where n is the number of simulations %interval= tinv(0.975,n-1)*S/sqrt(n), S is the estimation of the standard deviation
ACC_bin_l_up=ACC_bin_l_m+ACC_bin_l_int; % upper interval
ACC_bin_l_down=ACC_bin_l_m-ACC_bin_l_int; %lower interval

%binary; approx. linear; CDT; binary
ACC_spr_l_m=mean(ACC_spr,3); % estimation of the mean
ACC_spr_l_std=std(ACC_spr,0,3); % estimation of the standard deviation
ACC_spr_l_int= tinv(0.975,n-1)*ACC_spr_l_std/sqrt(n); %In order to estimate confidence intervals with 95% one need a coefficient from student-t distribution and with n-1 degrees of freedom where n is the number of simulations %interval= tinv(0.975,n-1)*S/sqrt(n), S is the estimation of the standard deviation
ACC_spr_l_up=ACC_spr_l_m+ACC_spr_l_int; % upper interval
ACC_spr_l_down=ACC_spr_l_m-ACC_spr_l_int; %lower interval


compareHoloGN_nonlinear %does the experiments for approx. nonlinear
%binary; approx. nonlinear; majority rule; binary
ACC_bin_nl_m=mean(ACC_bin,3); % estimation of the mean
ACC_bin_nl_std=std(ACC_bin,0,3); % estimation of the standard deviation
ACC_bin_nl_int= tinv(0.975,n-1)*ACC_bin_nl_std/sqrt(n); %In order to estimate confidence intervals with 95% one need a coefficient from student-t distribution and with n-1 degrees of freedom where n is the number of simulations %interval= tinv(0.975,n-1)*S/sqrt(n), S is the estimation of the standard deviation
ACC_bin_nl_up=ACC_bin_nl_m+ACC_bin_nl_int; % upper interval
ACC_bin_nl_down=ACC_bin_nl_m-ACC_bin_nl_int; %lower interval

%binary; approx. nonlinear; CDT; binary
ACC_spr_nl_m=mean(ACC_spr,3); % estimation of the mean
ACC_spr_nl_std=std(ACC_spr,0,3); % estimation of the standard deviation
ACC_spr_nl_int= tinv(0.975,n-1)*ACC_spr_nl_std/sqrt(n); %In order to estimate confidence intervals with 95% one need a coefficient from student-t distribution and with n-1 degrees of freedom where n is the number of simulations %interval= tinv(0.975,n-1)*S/sqrt(n), S is the estimation of the standard deviation
ACC_spr_nl_up=ACC_spr_nl_m+ACC_spr_nl_int; % upper interval
ACC_spr_nl_down=ACC_spr_nl_m-ACC_spr_nl_int; %lower interval

%binary; orthogonal; majority rule; binary
compareHoloGN_dense_HoloGN %does the experiments
ACC_bin_hgn_m=mean(ACC_bin_hgn,3); % estimation of the mean
ACC_bin_hgn_std=std(ACC_bin_hgn,0,3); % estimation of the standard deviation
ACC_bin_hgn_int= tinv(0.975,n-1)*ACC_bin_hgn_std/sqrt(n); %In order to estimate confidence intervals with 95% one need a coefficient from student-t distribution and with n-1 degrees of freedom where n is the number of simulations %interval= tinv(0.975,n-1)*S/sqrt(n), S is the estimation of the standard deviation
ACC_bin_hgn_up=ACC_bin_hgn_m+ACC_bin_hgn_int; % upper interval
ACC_bin_hgn_down=ACC_bin_hgn_m-ACC_bin_hgn_int; %lower interval

%binary; orthogonal; CDT; binary
compareHoloGN_sparse_HoloGN%does the experiments
ACC_spr_hgn_m=mean(ACC_spr_hgn,3); % estimation of the mean
ACC_spr_hgn_std=std(ACC_spr_hgn,0,3); % estimation of the standard deviation
ACC_spr_hgn_int= tinv(0.975,n-1)*ACC_spr_hgn_std/sqrt(n); %In order to estimate confidence intervals with 95% one need a coefficient from student-t distribution and with n-1 degrees of freedom where n is the number of simulations %interval= tinv(0.975,n-1)*S/sqrt(n), S is the estimation of the standard deviation
ACC_spr_hgn_up=ACC_spr_hgn_m+ACC_spr_hgn_int; % upper interval
ACC_spr_hgn_down=ACC_spr_hgn_m-ACC_spr_hgn_int; %lower interval


%Plot all results
figure()
for i=1:5

    
    subplot(2,3,i)
    hold on
    
    %plot confidence intervals
    
    hpatch=patch([DENS fliplr(DENS)], [ ACC_bip_down(i,:) fliplr(ACC_bip_up(i,:))], 0.50*[1 1 1]); %bipolar; linear; MAP; integer
    set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.5);    

    hpatch=patch([DENS fliplr(DENS)], [ ACC_bin_l_down(i,:) fliplr(ACC_bin_l_up(i,:))], 0.50*[1 1 1]); %binary; approx. linear; majority rule; binary
    set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.5);       
    
    hpatch=patch([DENS fliplr(DENS)], [ ACC_spr_l_down(i,:) fliplr(ACC_spr_l_up(i,:))], 0.50*[1 1 1]); %binary; approx. linear; CDT; binary
    set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    hpatch=patch([DENS fliplr(DENS)], [ ACC_bin_nl_down(i,:) fliplr(ACC_bin_nl_up(i,:))], 0.50*[1 1 1]); %binary; approx. nonlinear; majority rule; binary
    set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.5);       
    
    hpatch=patch([DENS fliplr(DENS)], [ ACC_spr_nl_down(i,:) fliplr(ACC_spr_nl_up(i,:))], 0.50*[1 1 1]); %binary; approx. nonlinear; CDT; binary
    set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    hpatch=patch([DENS fliplr(DENS)], [ ACC_bin_hgn_down(i,:) fliplr(ACC_bin_hgn_up(i,:))], 0.50*[1 1 1]); %binary; orthogonal; majority rule; binary
    set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.5);       
    
    hpatch=patch([DENS fliplr(DENS)], [ ACC_spr_hgn_down(i,:) fliplr(ACC_spr_hgn_up(i,:))], 0.50*[1 1 1]); %binary; orthogonal; CDT; binary
    set(hpatch, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    
    line_thick=2;
    %plot mean values
    l1=plot(DENS, ACC_bip_m(i,:),'-.g','Linewidth',line_thick); %bipolar; linear; MAP; integer
    l2=plot(DENS, ACC_bin_l_m(i,:),'b','Linewidth',line_thick); %binary; approx. linear; majority rule; binary
    l3=plot(DENS, ACC_spr_l_m(i,:),'--b','Linewidth',line_thick); %binary; approx. linear; CDT; binary
    l4=plot(DENS, ACC_bin_nl_m(i,:),'k','Linewidth',line_thick); %binary; approx. nonlinear; majority rule; binary
    l5=plot(DENS, ACC_spr_nl_m(i,:),'--k','Linewidth',line_thick); %binary; approx. nonlinear; CDT; binary    
    l6=plot(DENS, ACC_bin_hgn_m(i,:),'r','Linewidth',line_thick); %binary; orthogonal; majority rule; binary
    l7=plot(DENS, ACC_spr_hgn_m(i,:),'--r','Linewidth',line_thick); %binary; orthogonal; CDT; binary     
    

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
                legend([l1 l2 l3 l4 l5 l6 l7],{'bipolar; linear; MAP; integer','binary; approx. linear; majority rule; binary','binary; approx. linear; CDT; binary','binary; approx. nonlinear; majority rule; binary','binary; approx. nonlinear; CDT; binary','binary; orthogonal; majority rule; binary','binary; orthogonal; CDT; binary'}); 
                %legend('show')
                %legend('low-dim representation', 'dense, 10 000 elements', 'sparse, 100 000 elements', 'sparse, 10 000 elements', 'sparse, 2048 elements')
    end
    
end

