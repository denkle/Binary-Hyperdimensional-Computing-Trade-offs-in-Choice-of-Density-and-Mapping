%clear all
d=10000;%dimensionality
M1=100; %number of non-zero elements
%M=1; %sequnce length
K=1; 

Mrange=300; %number of stored elements
ACCm_ter=zeros(1,Mrange);

for M=1:Mrange
p1=1-(1-M1/d)^M; %superposition vector
p2=1-(1-p1)^K; %several permuted vectors
pc=p1*p2; %CDT vector
dp_hit=d*pc/(p1/(M1/d) ); %dotproduct with CDT
dp_rej=(M1)*pc;
ACCm_ter(1,M)=getaccuracy_sdr (dp_hit,dp_rej,D);   
end


