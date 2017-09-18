function RES = computeCDT_experimental_theoretical(DENS)

  D=10000; % dimensionality of vectors
  %DENS=[0.01]; % density of vectors
  N=16; % number of tokens in the bundle
  
  density=DENS(1);
  HD=round(rand(N,D)>(1-density));
  T=[0:10]; % number of CDT steps
  
  sim=6;
  RES=zeros(length(T),sim);
  for j=3:sim
  HD=round(rand(N,D)>(1-density));
      
  for i=1:length(T)
    t=T(i);
    RES(i,1)=t;
    RES(i,2)=computeCDT_pr(density,N,t); % analytical calculations
    RES(i,j)=sum(computeCDT( double(sum(HD)>0)  ,t))/D; % experimental calculations

  end
  
  end
  figure(1)
  hold on
  plot(RES(:,1),RES(:,2),'b','LineWidth',3) % analytical curve
  
  %RESm=mean(RES(:,3:sim),2);
  %plot(RES(:,1),RESm,'--r','LineWidth',1 )
  for j=3:sim
    plot(RES(:,1),RES(:,j),'--r','LineWidth',0.5 ) % simulation surves
  end
  
  %plot(RES(:,1),RES(:,3),'--r.','MarkerSize',16 )
  
  
  

  
  
  
end



function [pt] = computeCDT_pr(p1,N,t) % implements analytical calculations for density after CDT
 
    ps1=1-(1-p1)^N;
    
    if t==0
        pt=ps1;
    elseif t==1
        pt=ps1^2;    
    else
        pt=ps1^2;
        for i=2:t
           pt=pt+(ps1-pt)*ps1;
        end
    end

end



function [ngram] = computeCDT(ngram,K) % implements experimental calculations for density after CDT
% Performs Context-Dependent Thinning operation
% 
%   Input:
%       binary vector before thinning
%
%   Output:
%        resultant HD-vector   

%Dimensionality of HD-vectors
DIM=length(ngram(1,:));

coef_shift=round(0.8*DIM);

%Efficient implementation of disjunctive superposotion under the whole array
superposition=ngram;

%Implementation of fixed number of permutations for thinning
ngram=zeros(1,DIM); % resultant vector
%K=2;% number of permutations, i.e. thinning operations

%TH=zeros(K,DIM);

if K==0
    ngram=superposition;
else

for i=1:K
    shifted=circshift(superposition, [0 coef_shift+i]); % permuted vector   
    thinned=and(superposition,shifted); %thinned vector

    ngram=or(ngram,thinned); %update of the resultant vector    
    %disp([sum(thinned),sum(Res)]);
end
end

ngram=double(ngram);

end