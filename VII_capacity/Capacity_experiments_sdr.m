
%clear all
d=10000;%dimensionality
sh=round(0.8*d); %offset for shift in CDT operation

M1=100; %number of non-zero elements
%M=1; %sequnce length
K=1; %number of CDT runs

sim=100; %number of simulation
Mrange=300; %number of stored elements
ACCm=zeros(1,Mrange);
parfor M=1:Mrange
disp(M)
ACC=zeros(1,sim  );
for k=1:sim
    %disp(k)
    HD=round(rand(D,d)<=(M1/d)); %item memory
    %Encoding
        ls=M; %length of string
        HGN=zeros(1,d);    
        SEQ=zeros(1,ls);
        for i=1:ls
            %disp(i)
            SEQ(1,i)=round(1+(D-1)*rand(1,1)); %random index
            VEC=circshift(HD(SEQ(1,i),:),[0,i]);%permutation
            HGN=HGN+VEC;
            %Clipping                     
        end
        HGN=double(HGN>0);
        R=zeros(1,d);
        for i=1:K
        R=R+circshift(HGN,[0,sh+i]);
        end
        R=double(R>0); %thinning vector
        HGN=double(and(R,HGN)); %Thinned vector
        
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
ACCm(1,M)=mean(ACC);

end
%disp(mean(ACC))



