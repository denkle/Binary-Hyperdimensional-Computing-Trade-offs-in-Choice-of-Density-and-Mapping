function [Res] = CDT(array)
% Performs Context-Dependent Thinning operation
% 
%




%Number of HD-vectors to bundle
MAX=length(array(:,1));

%Dimensionality of HD-vectors
DIM=length(array(1,:));

coef_shift=round(0.8*DIM);

%Efficient implementation of disjunctive superposotion under the whole array
superposition=sum(array,1);
superposition=double(superposition>=1);


%Implementation of fixed number of permutations for thinning
Res=zeros(1,DIM); % resultant vector
K=1;% number of permutations

for i=1:K
    shifted=circshift(superposition, [0 coef_shift+i]);
    thinned=and(superposition,shifted);
    Res=or(Res,thinned); %update of the resultant vector    
    %disp([sum(thinned),sum(Res)]);
end
Res=double(Res);

end

