%This script outputs the figure corresponding to Fig. 10 in the article
figure()
hold on


%Analytical capacity
D=27; %dictionary size
Capacity_analytical; 
plot(ACCm_ter,'k','LineWidth',1.5)

D=100; %dictionary size
Capacity_analytical; 
plot(ACCm_ter,'b','LineWidth',1.5)

D=1000; %dictionary size
Capacity_analytical; 
plot(ACCm_ter,'r','LineWidth',1.5)


%Experimental capacity
D=27; %dictionary size
Capacity_experiments_sdr; 
plot(ACCm,'k--','LineWidth',1.5)

D=100; %dictionary size
Capacity_experiments_sdr; 
plot(ACCm,'b--','LineWidth',1.5)

D=1000; %dictionary size
Capacity_experiments_sdr; 
plot(ACCm,'r--','LineWidth',1.5)


xlabel('Number of stored elements')
ylabel('Probability of the correct retrieval')
legend('D=27', 'D=100', 'D=1000')
