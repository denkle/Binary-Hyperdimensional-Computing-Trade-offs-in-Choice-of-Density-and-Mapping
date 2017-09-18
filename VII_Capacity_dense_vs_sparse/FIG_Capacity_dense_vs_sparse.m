%This script outputs the figure corresponding to Fig. 11 in the article


%%
figure() % plot results of simulations
hold on

%% Capacity dense binomial with majority rule

Drange=[1000,2000,3000,4000,5000,6000, 10000, 20000, 30000]; %dimensionality range

d=2000;%dimensionality
sh=round(0.8*d); %offset for CDT operation
D=27; %dictionary size

%M=1; 

HYPERGEOM=0; % if generate from hypergeometric distribution

sim=1000; %number of simulations

for j=1:length(Drange)
    d=Drange(j); %dimensionality
    M1=0.5*d; %number of non-zero elements
    
disp(d)    
Mrange=round(0.03*d); %number of stored elements in the vector

ACCm=zeros(1,Mrange); %accuracy of retrieval
for M=round(0.027*d):Mrange
    %disp(M)
ACC=zeros(1,sim  );
for k=1:sim
    %disp(k)
    
    if HYPERGEOM==1
        HD=zeros(D,d,'single'); %item memory
        for i=1:D
            ind=randperm(d,M1);
            HD(i,ind)=1; %populate item memory
        end    
    else
        HD=round(rand(D,d)<=(M1/d)); %populate item memory
    end
    
    %Encoding
        ls=M; %length of string
        HGN=zeros(1,d);    
        SEQ=zeros(1,ls);
        for i=1:ls
            %disp(i)
            SEQ(1,i)=round(1+(D-1)*rand(1,1)); %random index
            VEC=circshift(HD(SEQ(1,i),:),[0,i]);%permutation
            HGN=HGN+VEC; %update representations with a new token from the dictionary
                   
        end
        thr=ls;
        if mod(ls,2)==0    
            HGN=HGN+round(rand(1,d));
            thr=thr+1;
        end
        HGN=single(HGN>thr/2); %majority rule
        
        
        %Decoding
        SEQd=zeros(1,ls);
        for i=1:ls
            VEC=circshift(HGN,[0,-i]);  %permutation back 
            Dist=zeros(1,D);
            for ham=1:D        
                Dist(1,ham)=sum(xor(HD(ham,:),VEC));               
            end         
            [~,SEQd(1,i)]=min(Dist);
        end
        acc=sum(SEQ==SEQd);
        ACC(1,k)=(acc/ls);
    
end

%disp(mean(ACC))
ACCm(1,M)=mean(ACC);

%find maximal value in Mrange with accuracy higher than 99%

if ACCm(1,M)>=0.99
    capacity_db(1,j)=M;
end

end

end


plot(Drange, capacity_db, 'k', 'Linewidth', 2) % plot results


%% Capacity sparse binomial with CDT, T=1

Drange=[1000,2000,3000,4000,5000,6000, 10000, 20000, 30000]; %dimensionality range

d=2000;%dimensionality
sh=round(0.8*d); %offset for CDT operation
D=27; %dictionary size
M1=85; %number of non-zero elements
%M=1; 


HYPERGEOM=0; % if generate from hypergeometric distribution

sim=1000; %number of simulations

for j=1:length(Drange)
    d=Drange(j); %dimensionality
disp(d)    
Mrange=round(0.0050*d); %number of stored elements in the vector

ACCm=zeros(1,Mrange); %accuracy of retrieval
for M=round(0.0038*d):Mrange
ACC=zeros(1,sim  );
for k=1:sim
    %disp(k)
    
    if HYPERGEOM==1
        HD=zeros(D,d,'single'); %item memory
        for i=1:D
            ind=randperm(d,M1);
            HD(i,ind)=1; %populate item memory
        end    
    else
        HD=round(rand(D,d)<=(M1/d)); %populate item memory
    end
    
    %Encoding
        ls=M; %length of string
        HGN=zeros(1,d);    
        SEQ=zeros(1,ls);
        for i=1:ls
            %disp(i)
            SEQ(1,i)=round(1+(D-1)*rand(1,1)); %random index
            VEC=circshift(HD(SEQ(1,i),:),[0,i]);%permutation
            HGN=HGN+VEC; %update representations with a new token from the dictionary
            %Clipping                     
        end
        HGN=single(HGN>0); %NO CDT
        %1 CDT iteration
        R=circshift(HGN,[0,sh]); %thinning vector
        HGN=single(and(R,HGN)); %Thinned vector       
        
        %Decoding
        SEQd=zeros(1,ls);
        for i=1:ls
            VEC=circshift(HGN,[0,-i]);  %permutation back 
            Dist=HD*VEC';
            [~,SEQd(1,i)]=max(Dist);
        end
        acc=sum(SEQ==SEQd);
        ACC(1,k)=(acc/ls);
    
end

%disp(mean(ACC))
ACCm(1,M)=mean(ACC);

%find maximal value in Mrange with accuracy higher than 99%

if ACCm(1,M)>=0.99
    capacity_sbcdt(1,j)=M;
end

end

end


plot(Drange, capacity_sbcdt, '-.b', 'Linewidth', 2) % plot results




%% Capacity sparse hypergeometric with CDT, T=1

Drange=[1000,2000,3000,4000,5000,6000, 10000, 20000, 30000]; %dimensionality range

d=2000;%dimensionality
sh=round(0.8*d); %offset for CDT operation
D=27; %dictionary size
M1=85; %number of non-zero elements
%M=1; 


HYPERGEOM=1; % if generate from hypergeometric distribution

sim=1000; %number of simulations

for j=1:length(Drange)
    d=Drange(j); %dimensionality
disp(d)    
Mrange=round(0.012*d); %number of stored elements in the vector

ACCm=zeros(1,Mrange); %accuracy of retrieval
for M=round(0.0095*d):Mrange
ACC=zeros(1,sim  );
for k=1:sim
    %disp(k)
    
    if HYPERGEOM==1
        HD=zeros(D,d,'single'); %item memory
        for i=1:D
            ind=randperm(d,M1);
            HD(i,ind)=1; %populate item memory
        end    
    else
        HD=round(rand(D,d)<=(M1/d)); %populate item memory
    end
    
    %Encoding
        ls=M; %length of string
        HGN=zeros(1,d);    
        SEQ=zeros(1,ls);
        for i=1:ls
            %disp(i)
            SEQ(1,i)=round(1+(D-1)*rand(1,1)); %random index
            VEC=circshift(HD(SEQ(1,i),:),[0,i]);%permutation
            HGN=HGN+VEC; %update representations with a new token from the dictionary
            %Clipping                     
        end
        HGN=single(HGN>0); %NO CDT
        %1 CDT iteration
        R=circshift(HGN,[0,sh]); %thinning vector
        HGN=single(and(R,HGN)); %Thinned vector       
        
        %Decoding
        SEQd=zeros(1,ls);
        for i=1:ls
            VEC=circshift(HGN,[0,-i]);  %permutation back 
            Dist=HD*VEC';
            [~,SEQd(1,i)]=max(Dist);
        end
        acc=sum(SEQ==SEQd);
        ACC(1,k)=(acc/ls);
    
end

%disp(mean(ACC))
ACCm(1,M)=mean(ACC);

%find maximal value in Mrange with accuracy higher than 99%

if ACCm(1,M)>=0.99
    capacity_shcdt(1,j)=M;
end

end

end

plot(Drange, capacity_shcdt, ':k', 'Linewidth', 2) % plot results


%% Capacity sparse hypergeometric without CDT

Drange=[1000,2000,3000,4000,5000,6000, 10000, 20000, 30000]; %dimensionality range

d=2000;%dimensionality
sh=round(0.8*d); %offset for CDT operation
D=27; %dictionary size
M1=85; %number of non-zero elements
%M=1; 


HYPERGEOM=1; % if generate from hypergeometric distribution

sim=1000; %number of simulations


for j=1:length(Drange)
    d=Drange(j); %dimensionality
disp(d)    
Mrange=round(0.03*d); %number of stored elements in the vector

ACCm=zeros(1,Mrange); %accuracy of retrieval
for M=round(0.027*d):Mrange
ACC=zeros(1,sim  );
for k=1:sim
    %disp(k)
    
    if HYPERGEOM==1
        HD=zeros(D,d,'single'); %item memory
        for i=1:D
            ind=randperm(d,M1);
            HD(i,ind)=1; %populate item memory
        end    
    else
        HD=round(rand(D,d)<=(M1/d)); %populate item memory
    end
    
    %Encoding
        ls=M; %length of string
        HGN=zeros(1,d);    
        SEQ=zeros(1,ls);
        for i=1:ls
            %disp(i)
            SEQ(1,i)=round(1+(D-1)*rand(1,1)); %random index
            VEC=circshift(HD(SEQ(1,i),:),[0,i]);%permutation
            HGN=HGN+VEC; %update representations with a new token from the dictionary
            %Clipping                     
        end
        HGN=single(HGN>0); %NO CDT
        
        %Decoding
        SEQd=zeros(1,ls);
        for i=1:ls
            VEC=circshift(HGN,[0,-i]);  %permutation back 
            Dist=HD*VEC';
            [~,SEQd(1,i)]=max(Dist);
        end
        acc=sum(SEQ==SEQd);
        ACC(1,k)=(acc/ls);
    
end

%disp(mean(ACC))
ACCm(1,M)=mean(ACC);

%find maximal value in Mrange with accuracy higher than 99%

if ACCm(1,M)>=0.99
    capacity_bloom(1,j)=M;
end


end

end

plot(Drange, capacity_bloom, '--r', 'Linewidth', 2) % plot results


%% legend and axes for the figure
xlabel('Dimensionality of HD vectors')
ylabel('Capacity of the bundle vector')
legend('dense HD vector', 'sparse, T=1, binomial', 'sparse, T=1, hypergeometric', 'sparse, T=0, hypergeometric')
box on
grid on




