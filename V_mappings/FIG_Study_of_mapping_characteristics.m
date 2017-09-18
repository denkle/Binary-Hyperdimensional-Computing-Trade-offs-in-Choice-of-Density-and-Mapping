%This script outputs the figure corresponding to Fig. 3 in the article

clear;
holoGN; % preload functions needed for sumulations
seed=1;

D=10000; % set dimensionality of HD vectors
global codingMode;
codingMode = 'dense_binary';
N=1;
MAXL=21;


DENS=[0.1:0.1:0.5]; 
%Create string for later use in the legend
for i=1:length(DENS)
s=num2str(DENS(i));
LEG(i,1:length(s))=s;
end



%% Ideal linear mapping
HD_dens=zeros(MAXL+1,length(DENS));
for ii=1:length(DENS)
density=DENS(ii);

iM = initItemMemories (D, MAXL, N,seed,density);
HD=zeros(MAXL+1,N);

for j=1:N
    density_act=sum(iM(strcat(int2str(0),'x',int2str(j))));
    for i=0:MAXL
        key = strcat(int2str(i),'x',int2str(j));
        
        HD(i+1,j)=sum(and(iM(key),iM(strcat(int2str(0),'x',int2str(j))) ));%/density_act ;
    end
end
HD=mean(HD,2);
HD=HD/(density*D);
HD_dens(:,ii)=HD; %Normalized Hamming distances for each density
end
%plot results for the current mapping scheme
figure
subplot(2,2,1)
plot(0:20,HD_dens(1:21,:),'Linewidth',2)
legend(LEG);
grid on
ylim([0,1])
xlabel('Signal level, mV')
ylabel('Dot product normalized by M')
title('Linear Mapping, ideal')


%% 

%% Approximate linear mapping
HD_dens=zeros(MAXL+1,length(DENS));
for ii=1:length(DENS)
density=DENS(ii);

iM = initItemMemoriesLinear (D, MAXL, N,seed,density);
HD=zeros(MAXL+1,N);

for j=1:N
    density_act=sum(iM(strcat(int2str(0),'x',int2str(j))));
    for i=0:MAXL
        key = strcat(int2str(i),'x',int2str(j));
        
        HD(i+1,j)=sum(and(iM(key),iM(strcat(int2str(0),'x',int2str(j))) ));%/density_act ;
    end
end
HD=mean(HD,2);
HD=HD/(density*D);
HD_dens(:,ii)=HD; %Normalized Hamming distances for each density
end
%plot results for the current mapping scheme
subplot(2,2,3)
plot(0:20,HD_dens(1:21,:),'Linewidth',2)
%legend(LEG);
grid on
ylim([0,1])

xlabel('Signal level, mV')
ylabel('Dot product normalized by M')
title('Linear Mapping, approximate')


%% Nonlinear mapping 
HD_dens=zeros(MAXL+1,length(DENS));
for ii=1:length(DENS)
density=DENS(ii);

iM = initItemMemoriesNonLinear (D, MAXL, N,seed,density);
HD=zeros(MAXL+1,N);

for j=1:N
    density_act=sum(iM(strcat(int2str(0),'x',int2str(j))));
    for i=0:MAXL
        key = strcat(int2str(i),'x',int2str(j));
        
        HD(i+1,j)=sum(and(iM(key),iM(strcat(int2str(0),'x',int2str(j))) ));%/density_act ;
    end
end
HD=mean(HD,2);
HD=HD/(density*D);
HD_dens(:,ii)=HD; %Normalized Hamming distances for each density
end
%plot results for the current mapping scheme
subplot(2,2,4)
plot(0:20,HD_dens(1:21,:),'Linewidth',2)
%legend(LEG);
grid on
ylim([0,1])

xlabel('Signal level, mV')
ylabel('Dot product normalized by M')
title('Nonlinear Mapping, approximate')


%% Orthogonal mapping 
HD_dens=zeros(MAXL+1,length(DENS));
for ii=1:length(DENS)
density=DENS(ii);

iM = initItemMemoriesHoloGN (D, MAXL, N,seed,density);
HD=zeros(MAXL+1,N);

for j=1:N
    density_act=sum(iM(strcat(int2str(0),'x',int2str(j))));
    for i=0:MAXL
        key = strcat(int2str(i),'x',int2str(j));
        
        HD(i+1,j)=sum(and(iM(key),iM(strcat(int2str(0),'x',int2str(j))) ));%/density_act ;
    end
end
HD=mean(HD,2);
HD=HD/(density*D); %Normalized Hamming distances for each density
HD_dens(:,ii)=HD;
end
%plot results for the current mapping scheme
subplot(2,2,2)
plot(0:20,HD_dens(1:21,:),'Linewidth',2)
%legend(LEG);
grid on
ylim([0,1])

xlabel('Signal level, mV')
ylabel('Dot product normalized by M')
title('Orthogonal Mapping')
